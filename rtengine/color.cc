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

#include "rtengine.h"
#include "color.h"
#include "iccmatrices.h"
#include "mytime.h"
#include "sleef.c"

using namespace std;

namespace rtengine {

    extern const Settings* settings;
	
    LUTf Color::cachef;
    LUTf Color::gamma2curve;

    LUTf Color::gammatab;
    LUTf Color::igammatab_srgb;
    LUTf Color::gammatab_srgb;
  //  LUTf Color::igammatab_709;
//	LUTf Color::gammatab_709;
    LUTf Color::igammatab_55;
	LUTf Color::gammatab_55;
    LUTf Color::igammatab_4;
	LUTf Color::gammatab_4;
	
    LUTf Color::igammatab_26_11;
    LUTf Color::gammatab_26_11;
    LUTf Color::igammatab_24_17;
    LUTf Color::gammatab_24_17a;

    // Wikipedia sRGB: Unlike most other RGB color spaces, the sRGB gamma cannot be expressed as a single numerical value.
    // The overall gamma is approximately 2.2, consisting of a linear (gamma 1.0) section near black, and a non-linear section elsewhere involving a 2.4 exponent
    // and a gamma (slope of log output versus log input) changing from 1.0 through about 2.3.
    const double Color::sRGBGamma = 2.2;
    const double Color::sRGBGammaCurve = 2.4;

    const double Color::eps_max=580.40756;  //(MAXVALF* 216.0f/24389.0);
    const double Color::kappa=24389.0/27.0;//903.29630;

    const float Color::D50x=0.9642f; //0.96422;
    const float Color::D50z=0.8249f; //0.82521;
    const double Color::u0=4.0*D50x/(D50x+15+3*D50z);
    const double Color::v0=9.0/(D50x+15+3*D50z);
    const double Color::epskap=8.0;
    /*
     * Munsell Lch correction
     * Copyright (c) 2011  Jacques Desmis <jdesmis@gmail.com>
    */
    // Munsell Lch LUTf : 195 LUT
    // about 70% data are corrected with significative corrections
    // almost all data are taken for BG, YR, G excepted a few extreme values with a slight correction
    // No LUTf for BG and Y : low corrections
    // Only between 5B and 5PB for L > 40 : under very low corrections for L < 40

    //give hue in function of L and C : Munsell  correction
    LUTf Color::_4P10, Color::_4P20, Color::_4P30, Color::_4P40, Color::_4P50, Color::_4P60;
    LUTf Color::_1P10, Color::_1P20, Color::_1P30, Color::_1P40, Color::_1P50, Color::_1P60;
    LUTf Color::_10PB10, Color::_10PB20, Color::_10PB30, Color::_10PB40, Color::_10PB50, Color::_10PB60;
    LUTf Color::_9PB10, Color::_9PB20, Color::_9PB30, Color::_9PB40, Color::_9PB50, Color::_9PB60, Color::_9PB70, Color::_9PB80;
    LUTf Color::_75PB10, Color::_75PB20, Color::_75PB30, Color::_75PB40, Color::_75PB50, Color::_75PB60, Color::_75PB70, Color::_75PB80;
    LUTf Color::_6PB10, Color::_6PB20, Color::_6PB30, Color::_6PB40, Color::_6PB50, Color::_6PB60, Color::_6PB70, Color::_6PB80;
    LUTf Color::_45PB10, Color::_45PB20, Color::_45PB30, Color::_45PB40, Color::_45PB50, Color::_45PB60, Color::_45PB70, Color::_45PB80;
    LUTf Color::_3PB10, Color::_3PB20, Color::_3PB30, Color::_3PB40, Color::_3PB50, Color::_3PB60, Color::_3PB70, Color::_3PB80;
    LUTf Color::_15PB10, Color::_15PB20, Color::_15PB30, Color::_15PB40, Color::_15PB50, Color::_15PB60, Color::_15PB70, Color::_15PB80;
    LUTf Color::_05PB40, Color::_05PB50, Color::_05PB60, Color::_05PB70, Color::_05PB80;
    LUTf Color::_10B40, Color::_10B50, Color::_10B60, Color::_10B70, Color::_10B80;
    LUTf Color::_9B40, Color::_9B50, Color::_9B60, Color::_9B70, Color::_9B80;
    LUTf Color::_7B40, Color::_7B50, Color::_7B60, Color::_7B70, Color::_7B80;
    LUTf Color::_5B40, Color::_5B50, Color::_5B60, Color::_5B70, Color::_5B80;
    LUTf Color::_10YR20, Color::_10YR30, Color::_10YR40, Color::_10YR50, Color::_10YR60, Color::_10YR70, Color::_10YR80, Color::_10YR90;
    LUTf Color::_85YR20, Color::_85YR30, Color::_85YR40, Color::_85YR50, Color::_85YR60, Color::_85YR70, Color::_85YR80, Color::_85YR90;
    LUTf Color::_7YR30, Color::_7YR40, Color::_7YR50, Color::_7YR60, Color::_7YR70, Color::_7YR80;
    LUTf Color::_55YR30, Color::_55YR40, Color::_55YR50, Color::_55YR60, Color::_55YR70, Color::_55YR80, Color::_55YR90;
    LUTf Color::_4YR30, Color::_4YR40, Color::_4YR50, Color::_4YR60, Color::_4YR70, Color::_4YR80;
    LUTf Color::_25YR30, Color::_25YR40, Color::_25YR50, Color::_25YR60, Color::_25YR70;
    LUTf Color::_10R30, Color::_10R40, Color::_10R50, Color::_10R60, Color::_10R70;
    LUTf Color::_9R30, Color::_9R40, Color::_9R50, Color::_9R60, Color::_9R70;
    LUTf Color::_7R30, Color::_7R40, Color::_7R50, Color::_7R60, Color::_7R70;
    LUTf Color::_5R10, Color::_5R20, Color::_5R30;
    LUTf Color::_25R10, Color::_25R20, Color::_25R30;
    LUTf Color::_10RP10, Color::_10RP20, Color::_10RP30;
    LUTf Color::_7G30, Color::_7G40, Color::_7G50, Color::_7G60, Color::_7G70, Color::_7G80;
    LUTf Color::_5G30, Color::_5G40, Color::_5G50, Color::_5G60, Color::_5G70, Color::_5G80;
    LUTf Color::_25G30, Color::_25G40, Color::_25G50, Color::_25G60, Color::_25G70, Color::_25G80;
    LUTf Color::_1G30, Color::_1G40, Color::_1G50, Color::_1G60, Color::_1G70, Color::_1G80;
    LUTf Color::_10GY30, Color::_10GY40, Color::_10GY50, Color::_10GY60, Color::_10GY70, Color::_10GY80;
    LUTf Color::_75GY30, Color::_75GY40, Color::_75GY50, Color::_75GY60, Color::_75GY70, Color::_75GY80;
    LUTf Color::_5GY30, Color::_5GY40, Color::_5GY50, Color::_5GY60, Color::_5GY70, Color::_5GY80;

#ifdef _DEBUG
    MunsellDebugInfo::MunsellDebugInfo() {
    	reinitValues();
    }
	void MunsellDebugInfo::reinitValues() {
    	maxdhue[0]=maxdhue[1]=maxdhue[2]=maxdhue[3]=0.0f;
    	maxdhuelum[0]=maxdhuelum[1]=maxdhuelum[2]=maxdhuelum[3]=0.0f;
    	depass=depassLum=0;
	}
#endif


    void Color::init () {

        int maxindex = 65536;
        cachef(maxindex,0/*LUT_CLIP_BELOW*/);

        gamma2curve(maxindex,0);

        for (int i=0; i<maxindex; i++) {
            if (i>eps_max) {
                cachef[i] = 327.68*( exp(1.0/3.0 * log((double)i / MAXVALF) ));
            }
            else {
                cachef[i] = 327.68*((kappa*i/MAXVALF+16.0)/116.0);
            }
        }

        for (int i=0; i<maxindex; i++) {
            gamma2curve[i] = (gamma2(i/65535.0) * 65535.0);
        }

        /*******************************************/

        gammatab(65536,0);
        igammatab_srgb(65536,0);
        gammatab_srgb(65536,0);
      //  igammatab_709(65536,0);
      //  gammatab_709(65536,0);
        igammatab_55(65536,0);
        gammatab_55(65536,0);
        igammatab_4(65536,0);
        gammatab_4(65536,0);
		
        igammatab_26_11(65536,0);
        gammatab_26_11(65536,0);
        igammatab_24_17(65536,0);
        gammatab_24_17a(65536,0);

        for (int i=0; i<65536; i++)
            gammatab_srgb[i] = (65535.0 * gamma2 (i/65535.0));
        for (int i=0; i<65536; i++)
            igammatab_srgb[i] = (65535.0 * igamma2 (i/65535.0));
        for (int i=0; i<65536; i++)
            gammatab[i] = (65535.0 * pow (i/65535.0, 0.454545));

 /*       for (int i=0; i<65536; i++)
            gammatab_709[i] = (65535.0 * gamma709 (i/65535.0));
        for (int i=0; i<65536; i++)
            igammatab_709[i] = (65535.0 * igamma709 (i/65535.0));
*/			
        for (int i=0; i<65536; i++)
            gammatab_55[i] = (65535.0 * gamma55 (i/65535.0));
        for (int i=0; i<65536; i++)
            igammatab_55[i] = (65535.0 * igamma55 (i/65535.0));
        for (int i=0; i<65536; i++)
            gammatab_4[i] = (65535.0 * gamma4 (i/65535.0));
        for (int i=0; i<65536; i++)
            igammatab_4[i] = (65535.0 * igamma4 (i/65535.0));

        for (int i=0; i<65536; i++)
            gammatab_26_11[i] = (65535.0 * gamma26_11 (i/65535.0));
        for (int i=0; i<65536; i++)
            igammatab_26_11[i] = (65535.0 * igamma26_11 (i/65535.0));
        for (int i=0; i<65536; i++) {
			float j = (float)i / 65535.0f;	
            gammatab_24_17a[i] = gamma24_17(j);
			}
        for (int i=0; i<65536; i++)
            igammatab_24_17[i] = (65535.0 * igamma24_17 (i/65535.0));
			
        /*FILE* f = fopen ("c.txt", "wt");
        for (int i=0; i<256; i++)
        fprintf (f, "%g %g\n", i/255.0, clower (i/255.0, 2.0, 1.0));
        fclose (f);*/

        initMunsell();

        linearGammaTRC = cmsBuildGamma(NULL, 1.0);
    }

    void Color::cleanup () {
        if (linearGammaTRC) cmsFreeToneCurve(linearGammaTRC);
    }

    void Color::rgb2hsl(float r, float g, float b, float &h, float &s, float &l) {

        double var_R = double(r) / 65535.0;
        double var_G = double(g) / 65535.0;
        double var_B = double(b) / 65535.0;

        double m = min(var_R,var_G,var_B);
        double M = max(var_R,var_G,var_B);
        double C = M - m;

        double l_ = (M+m)/2.;
        l = float(l_);

        if (C<0.00001 && C>-0.00001) {  // no fabs, slow!
            h = 0.f;
            s = 0.f;
        }
        else {
            double h_;
            if (l_ <= 0.5)
                s = float( (M-m) / (M+m) );
            else
                s = float( (M-m) / (2.0-M-m) );

            if      ( var_R == M ) h_ =      (var_G - var_B)/C;
            else if ( var_G == M ) h_ = 2. + (var_B - var_R)/C;
            else                   h_ = 4. + (var_R - var_G)/C;
            h = float(h_ /= 6.0);

            if ( h < 0.f )  h += 1.f;
            if ( h > 1.f )  h -= 1.f;
        }
    }

    double Color::hue2rgb(double p, double q, double t){
        if (t < 0.) t += 6.;
        else if( t > 6.) t -= 6.;

        if      (t < 1.) return p + (q - p) * t;
        else if (t < 3.) return q;
        else if (t < 4.) return p + (q - p) * (4. - t);
        else             return p;
    }

    void Color::hsl2rgb (float h, float s, float l, float &r, float &g, float &b) {

        if (s == 0)
            r = g = b = 65535.0f * l; //  achromatic
        else {
            double m2;
            double h_ = double(h);
            double s_ = double(s);
            double l_ = double(l);

            if (l <= 0.5f)
                m2 = l_ * (1.0 + s_);
            else {
                m2 = l_ + s_ - l_ * s_;
            }

            double m1 = 2.0 * l_ - m2;

            r = float(65535.0 * hue2rgb (m1, m2, h_ * 6.0 + 2.0));
            g = float(65535.0 * hue2rgb (m1, m2, h_ * 6.0));
            b = float(65535.0 * hue2rgb (m1, m2, h_ * 6.0 - 2.0));
        }
    }


    void Color::rgb2hsv(float r, float g, float b, float &h, float &s, float &v) {
        double var_R = r / 65535.0;
        double var_G = g / 65535.0;
        double var_B = b / 65535.0;

        double var_Min = min(var_R,var_G,var_B);
        double var_Max = max(var_R,var_G,var_B);
        double del_Max = var_Max - var_Min;
        v = var_Max;
        if (del_Max<0.00001 && del_Max>-0.00001) {  // no fabs, slow!
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
        int i = (int)h1;  // floor() is very slow, and h1 is always >0
        float f = h1 - i; // fractional part of h

        float p = v * ( 1 - s );
        float q = v * ( 1 - s * f );
        float t = v * ( 1 - s * ( 1 - f ) );

        float r1,g1,b1;

        if      (i==1)    {r1 = q;  g1 = v;  b1 = p;}
        else if (i==2)    {r1 = p;  g1 = v;  b1 = t;}
        else if (i==3)    {r1 = p;  g1 = q;  b1 = v;}
        else if (i==4)    {r1 = t;  g1 = p;  b1 = v;}
        else if (i==5)    {r1 = v;  g1 = p;  b1 = q;}
        else /*i==(0|6)*/ {r1 = v;  g1 = t;  b1 = p;}

        r = ((r1)*65535.0);
        g = ((g1)*65535.0);
        b = ((b1)*65535.0);
    }

    // Function copied for speed concerns
    // Not exactly the same as above ; this one return a result in the [0.0 ; 1.0] range
    void Color::hsv2rgb01 (float h, float s, float v, float &r, float &g, float &b) {
        float h1 = h*6; // sector 0 to 5
        int i = int(h1);
        float f = h1 - i; // fractional part of h

        float p = v * ( 1 - s );
        float q = v * ( 1 - s * f );
        float t = v * ( 1 - s * ( 1 - f ) );

        if      (i==1)    {r = q;  g = v;  b = p;}
        else if (i==2)    {r = p;  g = v;  b = t;}
        else if (i==3)    {r = p;  g = q;  b = v;}
        else if (i==4)    {r = t;  g = p;  b = v;}
        else if (i==5)    {r = v;  g = p;  b = q;}
        else /*(i==0|6)*/ {r = v;  g = t;  b = p;}
    }

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
    void Color::xyz2Prophoto (float x, float y, float z, float &r, float &g, float &b) {
        r = ((prophoto_xyz[0][0]*x + prophoto_xyz[0][1]*y + prophoto_xyz[0][2]*z)) ;
        g = ((prophoto_xyz[1][0]*x + prophoto_xyz[1][1]*y + prophoto_xyz[1][2]*z)) ;
        b = ((prophoto_xyz[2][0]*x + prophoto_xyz[2][1]*y + prophoto_xyz[2][2]*z)) ;
    }
    void Color::Prophotoxyz (float r, float g, float b, float &x, float &y, float &z) {
        x = ((xyz_prophoto[0][0]*r + xyz_prophoto[0][1]*g + xyz_prophoto[0][2]*b)) ;
        y = ((xyz_prophoto[1][0]*r + xyz_prophoto[1][1]*g + xyz_prophoto[1][2]*b)) ;
        z = ((xyz_prophoto[2][0]*r + xyz_prophoto[2][1]*g + xyz_prophoto[2][2]*b)) ;
    }

    void Color::rgbxyz (float r, float g, float b, float &x, float &y, float &z, double xyz_rgb[3][3]) {
        x = ((xyz_rgb[0][0]*r + xyz_rgb[0][1]*g + xyz_rgb[0][2]*b)) ;
        y = ((xyz_rgb[1][0]*r + xyz_rgb[1][1]*g + xyz_rgb[1][2]*b)) ;
        z = ((xyz_rgb[2][0]*r + xyz_rgb[2][1]*g + xyz_rgb[2][2]*b)) ;
    }
	
	
    void Color::xyz2rgb (float x, float y, float z, float &r, float &g, float &b, double rgb_xyz[3][3]) {
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

    // same for float
    void Color::xyz2rgb (float x, float y, float z, float &r, float &g, float &b, float rgb_xyz[3][3]) {
        r = ((rgb_xyz[0][0]*x + rgb_xyz[0][1]*y + rgb_xyz[0][2]*z)) ;
        g = ((rgb_xyz[1][0]*x + rgb_xyz[1][1]*y + rgb_xyz[1][2]*z)) ;
        b = ((rgb_xyz[2][0]*x + rgb_xyz[2][1]*y + rgb_xyz[2][2]*z)) ;
    }

    void Color::trcGammaBW (float &r, float &g, float &b, float gammabwr, float gammabwg, float gammabwb) {
        // correct gamma for black and white image : pseudo TRC curve of ICC profil
        b/=65535.0f;
        b= pow (max(b,0.0f), gammabwb);
        b *=65535.0f;
        r/=65535.0f;
        r= pow (max(r,0.0f), gammabwr);
        r *=65535.0f;
        g/=65535.0f;
        g= pow (max(g,0.0f), gammabwg);
        g *=65535.0f;
    }

    /** @brief Compute the B&W constants for the B&W processing and its tool's GUI
     *
     * @param setting BlackWhite::setting
     * @param setting BlackWhite::filter
     */
    void Color::computeBWMixerConstants (const Glib::ustring &setting, const Glib::ustring &filter,  const Glib::ustring &algo, float &mixerRed, float &mixerGreen,
                  float &mixerBlue, float mixerOrange, float mixerYellow, float mixerCyan, float mixerPurple, float mixerMagenta,
                  bool autoc, bool complement, float &kcorec, double &rrm, double &ggm, double &bbm)
    {
        float somm;
        float som = mixerRed+mixerGreen+mixerBlue;
        // rM = mixerRed, gM = mixerGreen, bM = mixerBlue !

        //presets
        if     (setting=="RGB-Abs" || setting=="ROYGCBPM-Abs")
            kcorec=som/100.f;

        if (!autoc) {
            //if     (setting=="RGB-Abs" || setting=="ROYGCBPM-Abs")  {} //Keep the RGB mixer values as is!
            //else if(setting=="RGB-Rel" || setting=="ROYGCBPM-Rel")  {} //Keep the RGB mixer values as is!
            if     (setting=="NormalContrast")    { mixerRed=43.f ; mixerGreen=33.f;  mixerBlue=30.f;  }
            else if(setting=="Panchromatic")      { mixerRed=33.3f; mixerGreen=33.3f; mixerBlue=33.3f; }
            else if(setting=="HyperPanchromatic") { mixerRed=41.f ; mixerGreen=25.f;  mixerBlue=34.f;  }
            else if(setting=="LowSensitivity")    { mixerRed=27.f ; mixerGreen=27.f;  mixerBlue=46.f;  }
            else if(setting=="HighSensitivity")   { mixerRed=30.f ; mixerGreen=28.f;  mixerBlue=42.f;  }
            else if(setting=="Orthochromatic")    { mixerRed=0.f  ; mixerGreen=42.f;  mixerBlue=58.f;  }
            else if(setting=="HighContrast")      { mixerRed=40.f ; mixerGreen=34.f;  mixerBlue=60.f;  }
            else if(setting=="Luminance")         { mixerRed=30.f ; mixerGreen=59.f;  mixerBlue=11.f;  }
            else if(setting=="Landscape")         { mixerRed=66.f ; mixerGreen=24.f;  mixerBlue=10.f;  }
            else if(setting=="Portrait")          { mixerRed=54.f ; mixerGreen=44.f;  mixerBlue=12.f;  }
            else if(setting=="InfraRed")          { mixerRed=-40.f; mixerGreen=200.f; mixerBlue=-17.f; }
        }

        rrm=mixerRed;
        ggm=mixerGreen;
        bbm=mixerBlue;

        somm=mixerRed+mixerGreen+mixerBlue;
        mixerRed=mixerRed/somm; mixerGreen=mixerGreen/somm; mixerBlue=mixerBlue/somm;
        float koymcp=0.f;

        if(setting=="ROYGCBPM-Abs" || setting=="ROYGCBPM-Rel") {
            float obM=0.f;
            float ogM=0.f;
            float orM=0.f;

            float ybM=0.f;
            float yrM=0.f;
            float ygM=0.f;

            float mgM=0.f;
            float mrM=0.f;
            float mbM=0.f;

            float pgM=0.f;
            float prM=0.f;
            float pbM=0.f;

            float crM=0.f;
            float cgM=0.f;
            float cbM=0.f;
			//printf("mixred=%f\n",mixerRed);

            float fcompl = 1.f;
            if(complement && algo=="SP") fcompl = 3.f;//special
			else if(complement && algo=="LI") fcompl = 1.5f;//linear
            // ponderate filters: report to R=G=B=33
            // I ponder RGB channel, not only orange or yellow or cyan, etc...it's my choice !
            if(mixerOrange != 33) {
				if (algo=="SP") {//special
					if (mixerOrange >= 33) orM = fcompl*(mixerOrange*0.67f - 22.11f)/100.f; else orM = fcompl*(-0.3f*mixerOrange +9.9f)/100.f;
					if (mixerOrange >= 33) ogM = fcompl*(-0.164f*mixerOrange+5.412f)/100.f; else ogM = fcompl*(0.4f*mixerOrange-13.2f)/100.f;
			   }
			   else if (algo=="LI") {//linear
				orM = fcompl*(mixerOrange - 33.f)/100.f; 
				ogM = fcompl*(0.5f*mixerOrange-16.5f)/100.f;
				}
                if(complement) obM =(-0.492f*mixerOrange+16.236f)/100.f;
                mixerRed   += orM;
                mixerGreen += ogM;
                mixerBlue  += obM;
                koymcp += (orM+ogM+obM);
			//	printf("mixred+ORange=%f\n",mixerRed);
			
            }
            if(mixerYellow != 33) {
                if (algo=="SP") yrM = fcompl*(-0.134f*mixerYellow+4.422f)/100.f;//22.4
                else if (algo=="LI")yrM = fcompl*(0.5f*mixerYellow-16.5f)/100.f;//22.4
                ygM = fcompl*(0.5f  *mixerYellow-16.5f )/100.f;
                if(complement) ybM =(-0.492f*mixerYellow+16.236f)/100.f;
                mixerRed   += yrM;
                mixerGreen += ygM;
                mixerBlue  += ybM;
                koymcp += (yrM+ygM+ybM);
            }
            if(mixerMagenta != 33) {
			 if (algo=="SP"){
					if(mixerMagenta >= 33) mrM = fcompl*( 0.67f *mixerMagenta-22.11f)/100.f; else mrM = fcompl*(-0.3f*mixerMagenta +9.9f)/100.f;
					if(mixerMagenta >= 33) mbM = fcompl*(-0.164f*mixerMagenta+5.412f)/100.f; else mbM = fcompl*( 0.4f*mixerMagenta-13.2f)/100.f;
					}
			 else if (algo=="LI"){
				mrM = fcompl*(mixerMagenta-33.f)/100.f; 
				mbM = fcompl*(0.5f*mixerMagenta-16.5f)/100.f; 
				}
                if(complement) mgM =(-0.492f*mixerMagenta+16.236f)/100.f;
                mixerRed   += mrM;
                mixerGreen += mgM;
                mixerBlue  += mbM;
                koymcp += (mrM+mgM+mbM);
            }
            if(mixerPurple != 33) {
				if (algo=="SP") prM = fcompl*(-0.134f*mixerPurple+4.422f)/100.f;
                else if (algo=="LI")prM = fcompl*(0.5f*mixerPurple-16.5f)/100.f;
                pbM = fcompl*(0.5f*mixerPurple-16.5f)/100.f;
                if(complement) pgM = (-0.492f*mixerPurple+16.236f)/100.f;
                mixerRed   += prM;
                mixerGreen += pgM;
                mixerBlue  += pbM;
                koymcp += (prM+pgM+pbM);
            }
            if(mixerCyan != 33) {
				if (algo=="SP")cgM = fcompl*(-0.134f*mixerCyan +4.422f)/100.f;
                else if (algo=="LI")cgM = fcompl*(0.5f*mixerCyan -16.5f)/100.f;
                cbM = fcompl*(0.5f*mixerCyan-16.5f)/100.f;
                if(complement) crM = (-0.492f*mixerCyan+16.236f)/100.f;
                mixerRed   += crM;
                mixerGreen += cgM;
                mixerBlue  += cbM;
                koymcp += (crM+cgM+cbM);
            }
        }
        if(setting=="ROYGCBPM-Abs")
            kcorec = koymcp+som/100.f;
        //Color filters
        float filred,filgreen,filblue;
        filred=1.f;filgreen=1.f;filblue=1.f;
        if          (filter=="None")        {filred=1.f;  filgreen=1.f;  filblue=1.f;}
        else if     (filter=="Red")         {filred=1.f;  filgreen=0.05f;filblue=0.f;}
        else if     (filter=="Orange")      {filred=1.f;  filgreen=0.6f; filblue=0.f;}
        else if     (filter=="Yellow")      {filred=1.f;  filgreen=1.f;  filblue=0.05f;}
        else if     (filter=="YellowGreen") {filred=0.6f; filgreen=1.f;  filblue=0.3f;}
        else if     (filter=="Green")       {filred=0.2f; filgreen=1.f;  filblue=0.3f;}
        else if     (filter=="Cyan")        {filred=0.05f;filgreen=1.f;  filblue=1.f;}
        else if     (filter=="Blue")        {filred=0.f;  filgreen=0.05f;filblue=1.f;}
        else if     (filter=="Purple")      {filred=1.f;  filgreen=0.05f;filblue=1.f;}


        mixerRed   = mixerRed   * filred;
        mixerGreen = mixerGreen * filgreen;
        mixerBlue  = mixerBlue  * filblue;

        mixerRed   = mixerRed   / (mixerRed + mixerGreen + mixerBlue);
        mixerGreen = mixerGreen / (mixerRed + mixerGreen + mixerBlue);
        mixerBlue  = mixerBlue  / (mixerRed + mixerGreen + mixerBlue);
		if(filter!="None") {
			som = mixerRed+mixerGreen+mixerBlue;
			if(setting=="RGB-Abs" || setting=="ROYGCBPM-Abs") kcorec = kcorec*som;
			}		
		
    }

    void Color::calcGamma (double pwr, double ts, int mode, int imax, double &gamma0, double &gamma1, double &gamma2, double &gamma3, double &gamma4, double &gamma5) {
        //from Dcraw (D.Coffin)
        int i;
        double g[6], bnd[2]={0,0};

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
        float LL=L/327.68f;
        float aa=a/327.68f;
        float bb=b/327.68f;
        float fy = (0.00862069f * LL) + 0.137932f; // (L+16)/116
        float fx = (0.002f * aa) + fy;
        float fz = fy - (0.005f * bb);
        x = 65535.0f*f2xyz(fx)*D50x;
        z = 65535.0f*f2xyz(fz)*D50z;
        y=(LL>epskap) ? 65535.0f*fy*fy*fy : 65535.0f*LL/kappa;
    }

    void Color::XYZ2Lab(float X, float Y, float Z, float &L, float &a, float &b) {

        float x = X/D50x;
        float z = Z/D50z;
        float y= Y;
        float fx,fy,fz;

        fx = (x<65535.0f ? cachef[std::max(x,0.f)] : (327.68f*exp(log(x/MAXVALF)/3.0f )));
        fy = (y<65535.0f ? cachef[std::max(y,0.f)] : (327.68f*exp(log(y/MAXVALF)/3.0f )));
        fz = (z<65535.0f ? cachef[std::max(z,0.f)] : (327.68f*exp(log(z/MAXVALF)/3.0f )));

        L = (116.0f *  fy - 5242.88f); //5242.88=16.0*327.68;
        a = (500.0f * (fx - fy) );
        b = (200.0f * (fy - fz) );
    }

    void Color::Lab2Yuv(float L, float a, float b, float &Y, float &u, float &v) {
        float fy = (0.00862069 * L/327.68) + 0.137932; // (L+16)/116
        float fx = (0.002 * a/327.68) + fy;
        float fz = fy - (0.005 * b/327.68);
        float LL=L/327.68;

        float X = 65535.0*f2xyz(fx)*D50x;
       // Y = 65535.0*f2xyz(fy);
        float Z = 65535.0*f2xyz(fz)*D50z;
        Y=(LL/327.68f>epskap) ? 65535.0*fy*fy*fy : 65535.0*LL/kappa;

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

        float fx = (X<65535.0 ? cachef[X] : (327.68*exp(log(X/MAXVALF)/3.0 )));
        float fy = (Y<65535.0 ? cachef[Y] : (327.68*exp(log(Y/MAXVALF)/3.0 )));
        float fz = (Z<65535.0 ? cachef[Z] : (327.68*exp(log(Z/MAXVALF)/3.0 )));

        L = (116.0 * fy - 5242.88); //5242.88=16.0*327.68;
        a = (500.0 * (fx - fy) );
        b = (200.0 * (fy - fz) );
    }
	
    double Color::f2xyz(double f) {
        const double epsilonExpInv3 = 6.0/29.0;
        const double kappaInv = 27.0/24389.0;  // inverse of kappa

        return (f > epsilonExpInv3) ? f*f*f : (116. * f - 16.) * kappaInv;
	
    }

    /*
     * Gamut mapping algorithm
     * Copyright (c) 2010-2011  Emil Martinec <ejmartin@uchicago.edu>
     *
     * Solutions to scaling u and v to XYZ paralleliped boundaries
     * Some equations:
     *
     *    fu(X,Y,Z) = 4 X/(X + 15 Y + 3 Z);
     *    fv(X,Y,Z) = 9 Y/(X + 15 Y + 3 Z);
     *
     * take the plane spanned by X=a*Xr+b*Xg+c*Xb etc with one of a,b,c equal to 0 or 1,
     * and itersect with the line u0+lam*u, or in other words solve
     *
     *    u0+lam*u=fu(X,Y,Z)
     *    v0+lam*v=fv(X,Y,Z)
     *
     * The value of lam is the scale factor that takes the color to the gamut boundary
     * columns of the matrix p=xyz_rgb are RGB tristimulus primaries in XYZ
     * c is the color fixed on the boundary; and m=0 for c=0, m=1 for c=255
     */
    void Color::gamutmap(float &X, float &Y, float &Z, const double p[3][3]) 
    {
        float u = 4*X/(X+15*Y+3*Z)-u0;
        float v = 9*Y/(X+15*Y+3*Z)-v0;

        float lam[3][2];
        float lam_min = 1.0;

        for (int c=0; c<3; c++)
            for (int m=0; m<2; m++) {

                int c1=(c+1)%3;
                int c2=(c+2)%3;

                lam[c][m] = (-(p[0][c1]*p[1][c]*((-12 + 3*u0 + 20*v0)*Y + 4*m*65535*v0*p[2][c2])) +
                    p[0][c]*p[1][c1]*((-12 + 3*u0 + 20*v0)*Y + 4*m*65535*v0*p[2][c2]) -
                    4*v0*p[0][c1]*(Y - m*65535*p[1][c2])*p[2][c] + 4*v0*p[0][c]*(Y - m*65535*p[1][c2])*p[2][c1] -
                    (4*m*65535*v0*p[0][c2] - 9*u0*Y)*(p[1][c1]*p[2][c] - p[1][c]*p[2][c1]));

                lam[c][m] /= (3*u*Y*(p[0][c1]*p[1][c] - p[1][c1]*(p[0][c] + 3*p[2][c]) + 3*p[1][c]*p[2][c1]) +
                    4*v*(p[0][c1]*(5*Y*p[1][c] + m*65535*p[1][c]*p[2][c2] + Y*p[2][c] - m*65535*p[1][c2]*p[2][c]) -
                    p[0][c]*(5*Y*p[1][c1] + m*65535*p[1][c1]*p[2][c2] + Y*p[2][c1] - m*65535*p[1][c2]*p[2][c1]) +
                    m*65535*p[0][c2]*(p[1][c1]*p[2][c] - p[1][c]*p[2][c1])));

                if (lam[c][m]<lam_min && lam[c][m]>0) lam_min=lam[c][m];

            }

            u = u*lam_min + u0;
            v = v*lam_min + v0;

            X = (9*u*Y)/(4*v); 
            Z = (12 - 3*u - 20*v)*Y/(4*v);
    }

	void Color::skinred ( double J, double h, double sres, double Sp, float dred, float protect_red, int sk, float rstprotection, float ko, double &s)
	{
		float factorskin, factorsat,factor, factorskinext, interm;
		float scale = 100.0f/100.1f;//reduction in normal zone
		float scaleext=1.0f;//reduction in transition zone
		float protect_redh;
		float deltaHH=0.3f;//HH value transition : I have choice 0.3 radians
		float HH;
		bool doskin=false;
		//rough correspondence between h (JC) and H (lab) that has relatively little importance because transitions that blur the correspondence is not linear
		if     ((float)h>8.6f  && (float)h<=74.f ) {HH=(1.15f/65.4f)*(float)h-0.0012f;  doskin=true;}//H > 0.15   H<1.3
		else if((float)h>0.f   && (float)h<=8.6f ) {HH=(0.19f/8.6f )*(float)h-0.04f;    doskin=true;}//H>-0.04 H < 0.15
		else if((float)h>355.f && (float)h<=360.f) {HH=(0.11f/5.0f )*(float)h-7.96f;    doskin=true;}//H>-0.15 <-0.04
		else if((float)h>74.f  && (float)h<95.f  ) {HH=(0.30f/21.0f)*(float)h+0.24285f; doskin=true;}//H>1.3  H<1.6

		if(doskin)
		{
			float chromapro=sres/Sp;
			if(sk==1){//in C mode to adapt dred to J
			if     (J<16.0)   dred = 40.0f;
			else if(J<22.0)   dred = (4.1666f)*(float)J -26.6f;
			else if(J<60.0)   dred = 55.0f;
			else if(J<70.0)   dred = -1.5f*(float)J +145.0f;
			else              dred = 40.0f;
			}
			if(chromapro>0.0) Color::scalered ( rstprotection, chromapro, 0.0, HH, deltaHH, scale, scaleext);//Scale factor
			if(chromapro>1.0) {interm=(chromapro-1.0f)*100.0f;
				factorskin= 1.0f+(interm*scale)/100.0f;
				factorskinext=1.0f+(interm*scaleext)/100.0f;}
			else {
				factorskin= chromapro ;
				factorskinext= chromapro ;
			}
			factorsat=chromapro;
			factor=factorsat;
			Color::transitred ( HH, s, dred, factorskin, protect_red, factorskinext, deltaHH, factorsat, factor);	//transition
			s*=factor;
		}
		else s=ko*sres;

	}
	void Color::skinredfloat ( float J, float h, float sres, float Sp, float dred, float protect_red, int sk, float rstprotection, float ko, float &s)
	{
		float factorskin, factorsat,factor, factorskinext, interm;
		float scale = 100.0f/100.1f;//reduction in normal zone
		float scaleext=1.0f;//reduction in transition zone
		float protect_redh;
		float deltaHH=0.3f;//HH value transition : I have choice 0.3 radians
		float HH;
		bool doskin=false;
		//rough correspondence between h (JC) and H (lab) that has relatively little importance because transitions that blur the correspondence is not linear
		if     ((float)h>8.6f  && (float)h<=74.f ) {HH=(1.15f/65.4f)*(float)h-0.0012f;  doskin=true;}//H > 0.15   H<1.3
		else if((float)h>0.f   && (float)h<=8.6f ) {HH=(0.19f/8.6f )*(float)h-0.04f;    doskin=true;}//H>-0.04 H < 0.15
		else if((float)h>355.f && (float)h<=360.f) {HH=(0.11f/5.0f )*(float)h-7.96f;    doskin=true;}//H>-0.15 <-0.04
		else if((float)h>74.f  && (float)h<95.f  ) {HH=(0.30f/21.0f)*(float)h+0.24285f; doskin=true;}//H>1.3  H<1.6

		if(doskin)
		{
			float chromapro=sres/Sp;
			if(sk==1){//in C mode to adapt dred to J
			if     (J<16.0)   dred = 40.0f;
			else if(J<22.0)   dred = (4.1666f)*(float)J -26.6f;
			else if(J<60.0)   dred = 55.0f;
			else if(J<70.0)   dred = -1.5f*(float)J +145.0f;
			else              dred = 40.0f;
			}
			if(chromapro>0.0) Color::scalered ( rstprotection, chromapro, 0.0, HH, deltaHH, scale, scaleext);//Scale factor
			if(chromapro>1.0) {interm=(chromapro-1.0f)*100.0f;
				factorskin= 1.0f+(interm*scale)/100.0f;
				factorskinext=1.0f+(interm*scaleext)/100.0f;}
			else {
				factorskin= chromapro ;
				factorskinext= chromapro ;
			}
			factorsat=chromapro;
			factor=factorsat;
			Color::transitred ( HH, s, dred, factorskin, protect_red, factorskinext, deltaHH, factorsat, factor);	//transition
			s*=factor;
		}
		else s=ko*sres;

	}



    void Color::SkinSatcdbl (float lum, float hue, float chrom, float skinprot, float &scale, bool ciec, bool neg, float b_l, float t_l, float t_r, float b_r, int choice) {
 
 
 
        float C9=0.0f, C8=0.0f, C7=0.0f, C4=0.0f, C3=0.0f, C2=0.0f, C1=0.0f;
        float H9=0.0f, H8=0.0f, H7=0.0f, H4=0.0f, H3=0.0f, H2=0.0f, H1=0.0f, H10=0.0f,H11=0.0f;
        H9=0.05f;H8=0.25f;H7=0.1f;H4=0.02f;H3=0.02f;H2=0.1f;H1=0.1f;H10=-0.2f;H11=-0.2f;//H10 and H11 are curious...H11=-0.8 ??
        C9=8.0f;C8=15.0f;C7=12.0f;C4=7.0f;C3=5.0f;C2=5.0f;C1=5.0f;

		if (ciec) {
		float HH;
		bool doskin;
		if     ((float)hue>8.6f  && (float)hue<=74.f ) {HH=(1.15f/65.4f)*(float)hue-0.0012f;  doskin=true;}//H > 0.15   H<1.3
		else if((float)hue>0.f   && (float)hue<=8.6f ) {HH=(0.19f/8.6f )*(float)hue-0.04f;    doskin=true;}//H>-0.04 H < 0.15
		else if((float)hue>355.f && (float)hue<=360.f) {HH=(0.11f/5.0f )*(float)hue-7.96f;    doskin=true;}//H>-0.15 <-0.04
		else if((float)hue>74.f  && (float)hue<95.f  ) {HH=(0.30f/21.0f)*(float)hue+0.24285f; doskin=true;}//H>1.3  H<1.6
		else if((float)hue>=95.f && (float)hue<137.5f) {HH= 0.01882*(float)hue-0.18823;}// H>1.6 H<2.4
		else if((float)hue>285.f && (float)hue<=355.f)  {HH=0.1642*(float)hue -5.982;}//HH>-1.3  HH <-0.15
		
		hue=HH;
		
		}
        // wide area for transition
		if((t_r-t_l)<0.55f) t_l=t_r+0.55f;//avoid too small range

        if      (lum >= 92.0f && (hue > b_l && hue < t_r) && (chrom > 7.0f && chrom < (18.0f))) scale = (100.f-skinprot*0.4f)/100.1f;
        else if (lum >= 85.0f && lum < 92.0f && (hue > b_l+0.05f && hue < t_r) && (chrom > 7.0f && chrom < (35.0f+C9))) scale = (100.f-skinprot*0.4f)/100.1f;
        else if ((lum >= 20.f && lum < 85.f) && (hue > (b_l+0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f+C9) )) scale = (100.f-skinprot*0.4f)/100.1f;
        else if (lum < 20.0f && (hue > (b_l+0.07f+H11) && hue < t_r-0.1f) && (chrom > 7.0f && chrom < (45.0f+C1) )) scale = (100.f-skinprot*0.4f)/100.1f;

        // wide area  skin color, useful if not accurate colorimetry or if the user has changed hue and saturation

        if      (lum >= 92.0f  && (hue > t_l+0.4f && hue < t_r) && (chrom > 7.0f && chrom < (15.0f))) scale = (100.f-skinprot*0.6f)/100.1f;
        else if (lum >= 85.0f && lum < 92.0f  && (hue > t_l+0.4f && hue < t_r-0.3f) && (chrom > 7.0f && chrom < (26.0f+C9))) scale = (100.f-skinprot*0.6f)/100.1f;
        else if ((lum >= 20.f && lum < 85.f) && (hue > (b_l+0.07f + H11) && hue < t_r-0.2f) && (chrom > 7.0f && chrom < (48.0f+C9) )) scale = (100.f-skinprot*0.6f)/100.1f;
        else if (lum < 20.0f  && (hue > (b_l+0.07f+H11) && hue < t_r-0.2f) && (chrom > 7.0f && chrom < (35.0f+C1) )) scale = (100.f-skinprot*0.6f)/100.1f;

        // "real" skin color : take into account a slightly usage of contrast and saturation in RT if option "skin" = 1

        if       (lum >= 85.0f  && (hue > (t_l+0.53f-H9) && hue < (t_r+H9)) && (chrom > 8.0f && chrom < (14.0f+C9))) scale = (100.f-skinprot)/100.1f;
        else if ((lum >= 70.0f && lum < 85.0f)  && (hue > t_l+0.15f && hue < (t_r-0.2f+H8)) && (chrom > 8.0f && chrom < (35.0f+C8))) scale = (100.f-skinprot)/100.1f;
        else if ((lum >= 52.0f && lum < 70.0f)  && (hue > t_l && hue < (t_r+H7)) && (chrom > 11.0f && chrom < (35.0f+C7))) scale = (100.f-skinprot)/100.1f;
        else if ((lum >= 35.0f && lum < 52.0f)  && (hue > t_l && hue <  (t_r+H4)) && (chrom > 13.0f && chrom < (37.0f+C4))) scale = (100.f-skinprot)/100.1f;
        else if ((lum >= 20.0f && lum < 35.0f)  && (hue > t_l && hue < (t_r+H3)) && (chrom > 7.0f && chrom <(35.0f+C3) )) scale = (100.f-skinprot)/100.1f;
        else if ((lum > 10.0f && lum < 20.0f)  && (hue > (t_l-0.25f + H10) && hue < (t_r-0.3f +H2)) && (chrom > 8.0f && chrom < (23.0f+C2))) scale = (100.f-skinprot)/100.1f;
        else if ((lum < 10.0f)  && (hue > (t_l -0.2f + H10) && hue < (t_r-0.3f+H1)) && (chrom > 8.0f && chrom < (23.0f+C1))) scale = (100.f-skinprot)/100.1f;
		
		//extended zone for hair, beard and if user adjust high value for skinprot
		if(skinprot > 85.f && chrom < 20.f && neg) {
			float modula = -0.0666f*skinprot + 6.66f;
			scale *= modula;
		}
    }

	

	void Color::scalered ( float rstprotection, float param, float limit, float HH, float deltaHH, float &scale,float &scaleext)
	{
		if(rstprotection<99.9999) {
			if(param > limit)
				scale = rstprotection/100.1f;
			if((HH< (1.3f+deltaHH) && HH >=1.3f))
				scaleext=HH*(1.0f-scale)/deltaHH + 1.0f - (1.3f+deltaHH)*(1.0f-scale)/deltaHH;    //transition for Hue (red - yellow)
			else if((HH< 0.15f && HH >(0.15f-deltaHH)))
				scaleext=HH*(scale-1.0f)/deltaHH + 1.0f - (0.15f-deltaHH)*(scale-1.0f)/deltaHH;   //transition for hue (red purple)
		}
	}
	
	void Color::transitred (float HH, float Chprov1, float dred, float factorskin, float protect_red, float factorskinext, float deltaHH, float factorsat, float &factor)
	{
		if(HH>=0.15f && HH<1.3f) {
			if (Chprov1<dred)
				factor = factorskin;
			else if(Chprov1<(dred+protect_red))
				factor = (factorsat-factorskin)/protect_red*Chprov1+factorsat-(dred+protect_red)*(factorsat-factorskin)/protect_red;
		}
		// then test if chroma is in the extanded range
		else if ( HH>(0.15f-deltaHH) || HH<(1.3f+deltaHH) ) {
			if (Chprov1 < dred)
				factor = factorskinext;// C=dred=55 => real max of skin tones
			else if (Chprov1 < (dred+protect_red))// transition
				factor = (factorsat-factorskinext)/protect_red*Chprov1+factorsat-(dred+protect_red)*(factorsat-factorskinext)/protect_red;
		}
	}
	
    /*
     * AllMunsellLch correction
     * Copyright (c) 2012  Jacques Desmis <jdesmis@gmail.com>
     *
     * This function corrects the color (hue) for changes in chromaticity and luminance
     * to use in a "for" or "do while" statement
     *
     * Parameters:
     *    bool lumaMuns : true => luminance correction (for delta L > 10) and chroma correction ; false : only chroma
     *    float Lprov1 , Loldd : luminance after and before
     *    float HH: hue before
     *    float Chprov1, CC : chroma after and before
     *    float coorectionHuechroma : correction Hue for chromaticity (saturation)
     *    float correctlum : correction Hue for luminance (brigtness, contrast,...)
     *    MunsellDebugInfo* munsDbgInfo: (Debug target only) object to collect information.
     */
#ifdef _DEBUG
    void Color::AllMunsellLch(bool lumaMuns, float Lprov1,float Loldd,float HH,float Chprov1,float CC,float &correctionHuechroma,float &correctlum, MunsellDebugInfo* munsDbgInfo)
#else
    void Color::AllMunsellLch(bool lumaMuns, float Lprov1, float Loldd, float HH, float Chprov1, float CC, float &correctionHuechroma, float &correctlum)
#endif
    {

        bool contin1,contin2;
        float correctionHue=0.0,correctionHueLum=0.0;
        float correctlumprov=0.0;
        float correctlumprov2=0.0;
        bool correctL=false;
        float huelimit[8]={-2.48,-0.55,0.44,1.52,1.87,3.09,-0.27,0.44};//limits hue of blue-purple, red-yellow, green-yellow, red-purple
        if(CC >= 6.0 && CC < 140) {          //if C > 140 we say C=140 (only in Prophoto ...with very large saturation)
            if (Chprov1 > 140) Chprov1=139;  //limits of LUTf
            if (Chprov1 < 6) Chprov1=6;
            for(int zo=1;zo<=4;zo++) {
                if(HH>huelimit[2*zo-2] && HH<huelimit[2*zo-1]) {
                    //zone=zo;
                    contin1=contin2=false;
                    correctL=false;
                    MunsellLch (Lprov1, HH,Chprov1, CC, correctionHue, zo, correctionHueLum, correctL);        //munsell chroma correction
#ifdef _DEBUG
                    float absCorrectionHue = fabs(correctionHue);
                    if(correctionHue !=0.0) {
                        int idx=zo-1;
#pragma omp critical (maxdhue)
{
munsDbgInfo->maxdhue[idx] = MAX(munsDbgInfo->maxdhue[idx], absCorrectionHue);
}
                    }
                    if(absCorrectionHue > 0.45)
#pragma omp atomic
                        munsDbgInfo->depass++;        //verify if no bug in calculation
#endif
                    correctionHuechroma=correctionHue;    //preserve
                    if(lumaMuns) {
                        if(correctL) {
                            //for Munsell luminance correction
                            correctlumprov=correctionHueLum;
                            contin1=true;
                            correctL=false;
                        }
                        correctionHueLum=0.0;
                        correctionHue=0.0;
                        if(fabs(Lprov1-Loldd) > 6.0) {
                            // correction if delta L significative..Munsell luminance
                            MunsellLch (Loldd, HH,Chprov1, Chprov1, correctionHue, zo, correctionHueLum, correctL);

                            if(correctL) {
                                correctlumprov2=correctionHueLum;
                                contin2=true;
                                correctL=false;
                            }
                            correctionHueLum=0.0;

                            if(contin1==true && contin2==true)
                                correctlum=correctlumprov2-correctlumprov;
#ifdef _DEBUG
                            float absCorrectLum = fabs(correctlum);
                            if(correctlum !=0.0) {
                                int idx=zo-1;
#pragma omp critical (maxdhuelum)
{
                                munsDbgInfo->maxdhuelum[idx] = MAX(munsDbgInfo->maxdhuelum[idx],absCorrectLum);
}
                            }
                            if(absCorrectLum > 0.35)
#pragma omp atomic
                                munsDbgInfo->depassLum++;    //verify if no bug in calculation
#endif
                        }
                    }
                }
            }

        }

#ifdef _DEBUG
        if     (correctlum < -0.35f) correctlum =-0.35f;
        else if(correctlum >  0.35f) correctlum = 0.35f;
        if     (correctionHuechroma<-0.45f) correctionHuechroma=-0.45f;
        else if(correctionHuechroma> 0.45f) correctionHuechroma= 0.45f;
#endif

    }

    /*
     * GamutLchonly correction
     * Copyright (c)2012  Jacques Desmis <jdesmis@gmail.com> and Jean-Christophe Frisch <natureh@free.fr>
     *
     * This function puts the data (Lab) in the gamut of "working profile":
     * it returns the corrected values of the chromaticity and luminance
     *
     * float HH : hue
     * float Lprov1 : input luminance value, sent back corrected
     * float Chprov1: input chroma value, sent back corrected
     * float R,G,B : red, green and blue value of the corrected color
     * double wip : working profile
     * bool isHLEnabled : if "highlight reconstruction " is enabled
     * float coef : a float number between [0.95 ; 1.0[... the nearest it is from 1.0, the more precise it will be... and the longer too as more iteration will be necessary)
     * bool neg and moreRGB : only in DEBUG mode to calculate iterations for negatives values and > 65535
     */
#ifdef _DEBUG
    void Color::gamutLchonly (float HH, float &Lprov1, float &Chprov1, float &R, float &G, float &B, double wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef, bool &neg, bool &more_rgb)
#else
    void Color::gamutLchonly (float HH, float &Lprov1, float &Chprov1, float &R, float &G, float &B, double wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef)
#endif
    {
        const float ClipLevel = 65535.0f;
        bool inGamut;
#ifdef _DEBUG
        neg=false, more_rgb=false;
#endif
        float2  sincosval = xsincosf(HH);
        do {
            inGamut=true;

            //Lprov1=LL;
            float aprov1=Chprov1*sincosval.y;
            float bprov1=Chprov1*sincosval.x;

            //conversion Lab RGB to limit Lab values - this conversion is useful before Munsell correction
            float fy = (0.00862069f *Lprov1 )+ 0.137932f;
            float fx = (0.002f * aprov1) + fy;
            float fz = fy - (0.005f * bprov1);

            float x_ = 65535.0f * f2xyz(fx)*D50x;
           // float y_ = 65535.0f * f2xyz(fy);
            float z_ = 65535.0f * f2xyz(fz)*D50z;
            float y_=(Lprov1>epskap) ? 65535.0*fy*fy*fy : 65535.0*Lprov1/kappa;

            xyz2rgb(x_,y_,z_,R,G,B,wip);

            // gamut control before saturation to put Lab values in future gamut, but not RGB
            if (R<0.0f || G<0.0f || B<0.0f) {
#ifdef _DEBUG
                neg=true;
#endif
                if (Lprov1 < 0.01f)
					Lprov1 = 0.01f;
                Chprov1 *= higherCoef; // decrease the chromaticity value
                if (Chprov1 <= 3.0f)
					Lprov1 += lowerCoef;
                inGamut = false;
			} else if (!isHLEnabled && (R>ClipLevel || G>ClipLevel || B>ClipLevel)) {

				// if "highlight reconstruction" is enabled or the point is completely white (clipped, no color), don't control Gamut
#ifdef _DEBUG
				more_rgb=true;
#endif
				if (Lprov1 > 99.999f)
					Lprov1 = 99.98f;
				Chprov1 *= higherCoef;
				if (Chprov1 <= 3.0f)
					Lprov1 -= lowerCoef;
				inGamut = false;
			}
		}
        while (!inGamut);
        //end first gamut control
    }


    /*
     * LabGamutMunsell
     * Copyright (c) 2012  Jacques Desmis <jdesmis@gmail.com>
     *
     * This function is the overall Munsell's corrections, but only on global statement: I think it's better to use local statement with AllMunsellLch
     * not for use in a "for" or "do while" loop
     * they are named accordingly :  gamutLchonly and AllMunsellLch
     * it can be used before and after treatment (saturation, gamma, luminance, ...)
     *
     * Parameters:
     *    Labimage *lab   :       RT Lab data
     *    float *Lold     :       luminance before   - data: i*width + j
     *    float *Cold     :       chrominance before - data: i*width + j
     *    bool corMuns    :       performs Munsell correction
     *    bool lumaMuns   :       (used only if corMuns=true)
     *                            true:  apply luma + chroma Munsell correction if delta L > 10;
     *                            false: only chroma correction only
     *    bool gamut            : performs gamutLch
     *    Glib::ustring &working: working profile's name
     *    bool multiThread      : parallelize the loop
     */
    void Color::LabGamutMunsell(LabImage *lab, float *Lold, float *Cold, bool corMunsell, bool lumaMuns, bool isHLEnabled, bool gamut, const Glib::ustring &working, bool multiThread ) {
#ifdef _DEBUG
        MyTime t1e,t2e;
        t1e.set();
        int negat=0, moreRGB=0;
        MunsellDebugInfo* MunsDebugInfo=NULL;
        if (corMunsell)
            MunsDebugInfo = new MunsellDebugInfo();

    #pragma omp parallel default(shared) firstprivate(MunsDebugInfo) reduction(+: negat, moreRGB) if (multiThread)
#else
    #pragma omp parallel default(shared) if (multiThread)
#endif
    {

        //unsigned int  N = lab->W*lab->H;
        int width = lab->W, height = lab->H;
        TMatrix wiprof = iccStore->workingSpaceInverseMatrix (working);
        double wip[3][3] = {
            {wiprof[0][0],wiprof[0][1],wiprof[0][2]},
            {wiprof[1][0],wiprof[1][1],wiprof[1][2]},
            {wiprof[2][0],wiprof[2][1],wiprof[2][2]}
        };
        float correctlum=0.0;
        float correctionHuechroma=0.0;

#pragma omp for schedule(dynamic, 10)
        for (int i=0; i<height; i++)
            for (int j=0; j<width; j++) {
                float HH=xatan2f(lab->b[i][j],lab->a[i][j]);
                float Chprov1=sqrt(SQR(lab->a[i][j]/327.68f) + SQR(lab->b[i][j]/327.68f));
                float Lprov1=lab->L[i][j]/327.68f;
                float Loldd, Coldd;
                if(gamut) {
#ifdef _DEBUG
                    bool neg, more_rgb;
#endif
                    float R, G, B;

                    //gamut control : Lab values are in gamut
#ifdef _DEBUG
                    gamutLchonly(HH, Lprov1, Chprov1, R, G, B, wip, isHLEnabled, 0.15f, 0.96f, neg, more_rgb);
#else
                    gamutLchonly(HH, Lprov1, Chprov1, R, G, B, wip, isHLEnabled, 0.15f, 0.96f);
#endif

#ifdef _DEBUG
                    if(neg) negat++;
                    if(more_rgb) moreRGB++;
#endif
                }

                Loldd = Lold[i*width + j];
                Coldd = Cold[i*width + j];
                lab->L[i][j] = Lprov1*327.68f;
                correctionHuechroma = 0.0;
                correctlum = 0.0;

                if(corMunsell)
#ifdef _DEBUG
                    AllMunsellLch(lumaMuns, Lprov1, Loldd, HH, Chprov1, Coldd, correctionHuechroma, correctlum, MunsDebugInfo);
#else
                    AllMunsellLch(lumaMuns, Lprov1, Loldd, HH, Chprov1, Coldd, correctionHuechroma, correctlum);
#endif

                HH+=correctlum;        //hue Munsell luminance correction

				correctlum = 0.0f;
				float2 sincosval = xsincosf(HH+correctionHuechroma);
				lab->a[i][j] = Chprov1*sincosval.y*327.68f;
				lab->b[i][j] = Chprov1*sincosval.x*327.68f;

            }
    } // end of parallelization

#ifdef _DEBUG
        t2e.set();
        if (settings->verbose) {
            printf("Color::LabGamutMunsell (correction performed in %d usec):\n", t2e.etime(t1e));
            printf("   Gamut              : G1negat=%iiter G165535=%iiter \n",negat,moreRGB);
            if (MunsDebugInfo) {
                printf("   Munsell chrominance: MaxBP=%1.2frad  MaxRY=%1.2frad  MaxGY=%1.2frad  MaxRP=%1.2frad  depass=%i\n", MunsDebugInfo->maxdhue[0],    MunsDebugInfo->maxdhue[1],    MunsDebugInfo->maxdhue[2],    MunsDebugInfo->maxdhue[3],    MunsDebugInfo->depass);
                printf("   Munsell luminance  : MaxBP=%1.2frad  MaxRY=%1.2frad  MaxGY=%1.2frad  MaxRP=%1.2frad  depass=%i\n", MunsDebugInfo->maxdhuelum[0] ,MunsDebugInfo->maxdhuelum[1], MunsDebugInfo->maxdhuelum[2], MunsDebugInfo->maxdhuelum[3], MunsDebugInfo->depassLum);
            }
            else {
           	printf("   Munsell correction wasn't requested\n");
            }
        }
        if (MunsDebugInfo)
        	delete MunsDebugInfo;
#endif

    }

    /*
     * MunsellLch correction
     * Copyright (c) 2012  Jacques Desmis <jdesmis@gmail.com>
     *
     * Find the right LUT and calculate the correction
     */
    void Color::MunsellLch (float lum, float hue, float chrom, float memChprov, float &correction, int zone, float &lbe, bool &correctL) {

        int x = int(memChprov);
        int y = int(chrom);

        //begin PB correction + sky
        if(zone==1) {
            if(lum > 5.0) {
				if(lum <15.0) {
					if( (hue >= (_15PB10[x] - 0.035)) && (hue < (_15PB10[x] + 0.052) && x<=45)) {if(y>49) y=49;correction =  _15PB10[y] - _15PB10[x] ;lbe=_15PB10[y];correctL=true;}
					else if (( hue>=( _3PB10[x] -0.052))  && (hue < (_45PB10[x] + _3PB10[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB10[y] - _3PB10[x]; lbe =_3PB10[y];correctL=true;}
					else if (( hue>=(_45PB10[x] + _3PB10[x])/2.0)  && (hue < (_45PB10[x] +0.052)) && x <= 85) {if(y>89) y=89;correction =  _45PB10[y] - _45PB10[x] ;lbe=_45PB10[y];correctL=true;}
					else if (( hue>=(_6PB10[x] -0.052)  && (hue < (_6PB10[x] + _75PB10[x])/2.0))) {correction =  _6PB10[y] - _6PB10[x] ;lbe=_6PB10[y];correctL=true;}
					else if (( hue>=(_6PB10[x] + _75PB10[x])/2.0)  && (hue < (_9PB10[x] + _75PB10[x])/2.0)) {correction =  _75PB10[y] - _75PB10[x] ;lbe=_75PB10[y];correctL=true;}
					else if (( hue>=(_9PB10[x] + _75PB10[x])/2.0)  && (hue < (_9PB10[x] + _10PB10[x])/2.0)) {correction =  _9PB10[y] - _9PB10[x] ; lbe=_9PB10[y];correctL=true;}
					else if (( hue>=(_10PB10[x] + _9PB10[x])/2.0)  && (hue < (_1P10[x] + _10PB10[x])/2.0)) {correction =  _10PB10[y] - _10PB10[x] ;lbe=_10PB10[y];correctL=true;}
					else if (( hue>=(_10PB10[x] + _1P10[x])/2.0)  && (hue < (_1P10[x] + _4P10[x])/2.0)) {correction =  _1P10[y] - _1P10[x];lbe=_1P10[y];correctL=true;}
					else if (( hue>=(_1P10[x] + _4P10[x])/2.0)  && (hue < (0.035 + _4P10[x])/2.0)) {correction =  _4P10[y] - _4P10[x] ;lbe=_4P10[y];correctL=true;}
				}
				else if (lum <25.0) {
					if( (hue >= (_15PB20[x] - 0.035)) && (hue < (_15PB20[x] + _3PB20[x])/2.0) && x<=85) {if(y>89) y=89;correction =  _15PB20[y] - _15PB20[x] ;lbe= _15PB20[y];correctL=true;}
					else if (( hue>=(_15PB20[x] + _3PB20[x])/2.0)  && (hue < (_45PB20[x] + _3PB20[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB20[y] - _3PB20[x] ;lbe= _3PB20[y];correctL=true;}
					else if (( hue>=(_45PB20[x] + _3PB20[x])/2.0)  && (hue < ( _45PB20[x] + 0.052)) && x <= 85) {if(y>89) y=89;correction =  _45PB20[y] - _45PB20[x] ;lbe=_45PB20[y];correctL=true;}
					else if (( hue>=(_45PB20[x] + 0.052))  && (hue < (_6PB20[x] + _75PB20[x])/2.0)) {correction =  _6PB20[y] - _6PB20[x];lbe=_6PB20[y];correctL=true;}
					else if (( hue>=(_6PB20[x] + _75PB20[x])/2.0)  && (hue < (_9PB20[x] + _75PB20[x])/2.0)) {correction =  _75PB20[y] - _75PB20[x] ;lbe=_75PB20[y];correctL=true;}
					else if (( hue>=(_9PB20[x] + _75PB20[x])/2.0)  && (hue < (_9PB20[x] + _10PB20[x])/2.0)) {correction =  _9PB20[y] - _9PB20[x] ;lbe= _9PB20[y];correctL=true; }
					else if (( hue>=(_10PB20[x] + _9PB20[x])/2.0)  && (hue < (_1P20[x] + _10PB20[x])/2.0)) {correction =  _10PB20[y] - _10PB20[x] ;lbe= _10PB20[y];correctL=true;}
					else if (( hue>=(_10PB20[x] + _1P20[x])/2.0)  && (hue < (_1P20[x] + _4P20[x])/2.0)) {correction =  _1P20[y] - _1P20[x] ; lbe=_1P20[y];correctL=true;}
					else if (( hue>=(_1P20[x] + _4P20[x])/2.0)  && (hue < (0.035 + _4P20[x])/2.0)) {correction =  _4P20[y] - _4P20[x] ; lbe=_4P20[y];correctL=true;}
				}
				else if (lum <35.0) {
					if( (hue >= (_15PB30[x] - 0.035)) && (hue < (_15PB30[x] + _3PB30[x])/2.0) && x<=85 ) {if(y>89) y=89;correction =  _15PB30[y] - _15PB30[x] ;lbe=_15PB30[y];correctL=true;}
					else if (( hue>=(_15PB30[x] + _3PB30[x])/2.0)  && (hue < (_45PB30[x] + _3PB30[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB30[y] - _3PB30[x] ;lbe=_3PB30[y];correctL=true;}
					else if (( hue>=(_45PB30[x] + _3PB30[x])/2.0)  && (hue < (_45PB30[x]+0.052)) && x <= 85) {if(y>89) y=89;correction =  _45PB30[y] - _45PB30[x] ;lbe= _45PB30[y];correctL=true;}
					else if (( hue>=( _45PB30[x]+ 0.052))  && (hue < (_6PB30[x] + _75PB30[x])/2.0)) {correction =  _6PB30[y] - _6PB30[x] ; lbe=_6PB30[y];correctL=true;}
					else if (( hue>=(_6PB30[x] + _75PB30[x])/2.0)  && (hue < (_9PB30[x] + _75PB30[x])/2.0)) {correction =  _75PB30[y] - _75PB30[x] ;lbe= _75PB30[y] ;correctL=true;}
					else if (( hue>=(_9PB30[x] + _75PB30[x])/2.0)  && (hue < (_9PB30[x] + _10PB30[x])/2.0)) {correction =  _9PB30[y] - _9PB30[x] ;lbe=_9PB30[y]; correctL=true;}
					else if (( hue>=(_10PB30[x] + _9PB30[x])/2.0)  && (hue < (_1P30[x] + _10PB30[x])/2.0)) {correction =  _10PB30[y] - _10PB30[x] ;lbe=_10PB30[y];correctL=true;}
					else if (( hue>=(_10PB30[x] + _1P30[x])/2.0)  && (hue < (_1P30[x] + _4P30[x])/2.0)) {correction =  _1P30[y] - _1P30[x] ;lbe=_1P30[y];correctL=true; }
					else if (( hue>=(_1P30[x] + _4P30[x])/2.0)  && (hue < (0.035 + _4P30[x])/2.0)) {correction =  _4P30[y] - _4P30[x] ;lbe=_4P30[y];correctL=true;}
				}
				else if (lum <45.0) {
					if( (hue <= (_05PB40[x] + _15PB40[x])/2.0) && (hue > (_05PB40[x] + _10B40[x])/2.0) && x<75 ) {if(y>75) y=75; correction =  _05PB40[y] - _05PB40[x] ;lbe=_05PB40[y];correctL=true;}
					else if( (hue <= (_05PB40[x] + _10B40[x])/2.0) && (hue >(_10B40[x] + _9B40[x])/2.0) && x<70 ) {if(y>70) y=70;correction =  _10B40[y] - _10B40[x] ;lbe=_10B40[y];correctL=true;}
					else if( (hue <= (_10B40[x] + _9B40[x])/2.0) && (hue >(_9B40[x] + _7B40[x])/2.0) && x<70 ) {if(y>70) y=70;correction =  _9B40[y] - _9B40[x] ;lbe=_9B40[y];correctL=true;}
					else if( (hue <= (_9B40[x] + _7B40[x])/2.0) && (hue >(_5B40[x] + _7B40[x])/2.0) && x<70 ) {if(y>70) y=70;correction =  _7B40[y] - _7B40[x] ;lbe= _7B40[y];correctL=true;}
					else if (( hue<=(_5B40[x] + _7B40[x])/2.0)  && (hue > (_5B40[x]-0.035)) && x < 70) {if(y>70) y=70; correction =  _5B40[y] - _5B40[x] ;lbe=  _5B40[y];correctL=true;}    //

					else if( (hue >= (_15PB40[x] - 0.035)) && (hue < (_15PB40[x] + _3PB40[x])/2.0) && x<=85 ) {if(y>89) y=89;correction =  _15PB40[y] - _15PB40[x] ; lbe= _15PB40[y];correctL=true;}
					else if (( hue>=(_15PB40[x] + _3PB40[x])/2.0)  && (hue < (_45PB40[x] + _3PB40[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB40[y] - _3PB40[x] ;lbe=_3PB40[y];correctL=true;}
					else if (( hue>=(_45PB40[x] + _3PB40[x])/2.0)  && (hue < (_45PB40[x]+0.052)) && x <= 85) {if(y>89) y=89;correction =  _45PB40[y] - _45PB40[x] ;lbe=_45PB40[y] ;correctL=true;}
					else if (( hue>=(_45PB40[x]+0.052))  && (hue < (_6PB40[x] + _75PB40[x])/2.0)) {correction =  _6PB40[y] - _6PB40[x] ;lbe=_6PB40[y];correctL=true; }
					else if (( hue>=(_6PB40[x] + _75PB40[x])/2.0)  && (hue < (_9PB40[x] + _75PB40[x])/2.0)) {correction =  _75PB40[y] - _75PB40[x] ; lbe=_75PB40[y];correctL=true;}
					else if (( hue>=(_9PB40[x] + _75PB40[x])/2.0)  && (hue < (_9PB40[x] + _10PB40[x])/2.0)) {correction =  _9PB40[y] - _9PB40[x] ;lbe= _9PB40[y]; correctL=true;}
					else if (( hue>=(_10PB40[x] + _9PB40[x])/2.0)  && (hue < (_1P40[x] + _10PB40[x])/2.0)) {correction =  _10PB40[y] - _10PB40[x] ;lbe=_10PB40[y];correctL=true;}
					else if (( hue>=(_10PB40[x] + _1P40[x])/2.0)  && (hue < (_1P40[x] + _4P40[x])/2.0)) {correction =  _1P40[y] - _1P40[x] ;lbe=_1P40[y];correctL=true;}
					else if (( hue>=(_1P40[x] + _4P40[x])/2.0)  && (hue < (0.035 + _4P40[x])/2.0)) {correction =  _4P40[y] - _4P40[x] ;lbe= _4P40[y];correctL=true;}
				}
				else if (lum <55.0) {
					if( (hue <= (_05PB50[x] + _15PB50[x])/2.0) && (hue > (_05PB50[x] + _10B50[x])/2.0) && x<79 ) {if(y>79) y=79; correction =  _05PB50[y] - _05PB50[x] ;lbe= _05PB50[y];correctL=true;}
					else if( (hue <= (_05PB50[x] + _10B50[x])/2.0) && (hue >(_10B50[x] + _9B50[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _10B50[y] - _10B50[x] ;lbe=_10B50[y];correctL=true;}
					else if( (hue <= (_10B50[x] + _9B50[x])/2.0) && (hue >(_9B50[x] + _7B50[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _9B50[y] - _9B50[x] ;lbe=_9B50[y];correctL=true;}
					else if( (hue <= (_9B50[x] + _7B50[x])/2.0) && (hue >(_5B50[x] + _7B50[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _7B50[y] - _7B50[x] ;lbe=_7B50[y];correctL=true;}
					else if (( hue<=(_5B50[x] + _7B50[x])/2.0)  && (hue > (_5B50[x]-0.035)) && x < 79) {if(y>79) y=79; correction =  _5B50[y] - _5B50[x] ;lbe=_5B50[y];correctL=true; }    //

					else if( (hue >= (_15PB50[x] - 0.035)) && (hue < (_15PB50[x] + _3PB50[x])/2.0) && x<=85 ) {if(y>89) y=89;correction =  _15PB50[y] - _15PB50[x] ; lbe= _15PB50[y];correctL=true;}
					else if (( hue>=(_15PB50[x] + _3PB50[x])/2.0)  && (hue < (_45PB50[x] + _3PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB50[y] - _3PB50[x] ;lbe=_3PB50[y];correctL=true;}
					else if (( hue>=(_45PB50[x] + _3PB50[x])/2.0)  && (hue < (_6PB50[x] + _45PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _45PB50[y] - _45PB50[x] ;lbe=_45PB50[y];correctL=true; }
					else if (( hue>=(_6PB50[x] + _45PB50[x])/2.0)  && (hue < (_6PB50[x] + _75PB50[x])/2.0) && x <=85) {if(y>89) y=89;correction =  _6PB50[y] - _6PB50[x] ;lbe=_6PB50[y];correctL=true;}
					else if (( hue>=(_6PB50[x] + _75PB50[x])/2.0)  && (hue < (_9PB50[x] + _75PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _75PB50[y] - _75PB50[x] ;lbe=_75PB50[y];correctL=true;}
					else if (( hue>=(_9PB50[x] + _75PB50[x])/2.0)  && (hue < (_9PB50[x] + _10PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _9PB50[y] - _9PB50[x] ;lbe=_9PB50[y];correctL=true;}
					else if (( hue>=(_10PB50[x] + _9PB50[x])/2.0)  && (hue < (_1P50[x] + _10PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _10PB50[y] - _10PB50[x] ;lbe=_10PB50[y];correctL=true;}
					else if (( hue>=(_10PB50[x] + _1P50[x])/2.0)  && (hue < (_1P50[x] + _4P50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _1P50[y] - _1P50[x] ;lbe=_1P50[y];correctL=true; }
					else if (( hue>=(_1P50[x] + _4P50[x])/2.0)  && (hue < (0.035 + _4P50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _4P50[y] - _4P50[x] ;lbe=_4P50[y];correctL=true;}
				}
				else if (lum <65.0) {
					if( (hue <= (_05PB60[x] + _15PB60[x])/2.0) && (hue > (_05PB60[x] + _10B60[x])/2.0) && x<79 ) {if(y>79) y=79; correction =  _05PB60[y] - _05PB60[x] ;lbe=_05PB60[y];correctL=true;}
					else if( (hue <= (_05PB60[x] + _10B60[x])/2.0) && (hue >(_10B60[x] + _9B60[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _10B60[y] - _10B60[x] ;lbe= _10B60[y];correctL=true;}
					else if( (hue <= (_10B60[x] + _9B60[x])/2.0) && (hue >(_9B60[x] + _7B60[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _9B60[y] - _9B60[x] ;lbe=_9B60[y];correctL=true;}
					else if( (hue <= (_9B60[x] + _7B60[x])/2.0) && (hue >(_5B60[x] + _7B60[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _7B60[y] - _7B60[x] ;lbe= _7B60[y];correctL=true;}
					else if (( hue<=(_5B60[x] + _7B60[x])/2.0)  && (hue > (_5B60[x]-0.035)) && x < 79) {if(y>79) y=79; correction =  _5B60[y] - _5B60[x] ;lbe= _5B60[y];correctL=true;}    //

					else if( (hue >= (_15PB60[x] - 0.035)) && (hue < (_15PB60[x] + _3PB60[x])/2.0) && x<=85 ) {if(y>89) y=89;correction =  _15PB60[y] - _15PB60[x] ;lbe=_15PB60[y];correctL=true; }
					else if (( hue>=(_15PB60[x] + _3PB60[x])/2.0)  && (hue < (_45PB60[x] + _3PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB60[y] - _3PB60[x] ;lbe=_3PB60[y];correctL=true;}
					else if (( hue>=(_45PB60[x] + _3PB60[x])/2.0)  && (hue < (_6PB60[x] + _45PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _45PB60[y] - _45PB60[x] ;lbe=_45PB60[y];correctL=true;}
					else if (( hue>=(_6PB60[x] + _45PB60[x])/2.0)  && (hue < (_6PB60[x] + _75PB60[x])/2.0) && x <=85) {if(y>89) y=89;correction =  _6PB60[y] - _6PB60[x] ;lbe= _6PB60[y];correctL=true;}
					else if (( hue>=(_6PB60[x] + _75PB60[x])/2.0)  && (hue < (_9PB60[x] + _75PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _75PB60[y] - _75PB60[x] ;lbe= _75PB60[y];correctL=true;}
					else if (( hue>=(_9PB60[x] + _75PB60[x])/2.0)  && (hue < (_9PB60[x] + _10PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _9PB60[y] - _9PB60[x] ;lbe= _9PB60[y];correctL=true;}
					else if (( hue>=(_10PB60[x] + _9PB60[x])/2.0)  && (hue < (_1P60[x] + _10PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _10PB60[y] - _10PB60[x] ;lbe=_10PB60[y]; correctL=true;}
					else if (( hue>=(_10PB60[x] + _1P60[x])/2.0)  && (hue < (_1P60[x] + _4P60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _1P60[y] - _1P60[x] ; lbe= _1P60[y];correctL=true;}
					else if (( hue>=(_1P60[x] + _4P60[x])/2.0)  && (hue < (0.035 + _4P60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _4P60[y] - _4P60[x] ;lbe=_4P60[y];correctL=true; }
				}
				else if (lum < 75.0) {
					if( (hue <= (_05PB70[x] + _15PB70[x])/2.0) && (hue > (_05PB70[x] + _10B70[x])/2.0) && x<50 ) {if(y>49) y=49; correction =  _05PB70[y] - _05PB70[x] ;lbe=_05PB70[y];correctL=true;}
					else if( (hue <= (_05PB70[x] + _10B70[x])/2.0) && (hue >(_10B70[x] + _9B70[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _10B70[y] - _10B70[x] ;lbe=_10B70[y];correctL=true;}
					else if( (hue <= (_10B70[x] + _9B70[x])/2.0) && (hue >(_9B70[x] + _7B70[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _9B70[y] - _9B70[x] ;lbe= _9B70[y];correctL=true;}
					else if( (hue <= (_9B70[x] + _7B70[x])/2.0) && (hue >(_5B70[x] + _7B70[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _7B70[y] - _7B70[x] ;lbe=_7B70[y];correctL=true;}
					else if (( hue<=(_5B70[x] + _7B70[x])/2.0)  && (hue > (_5B70[x]-0.035)) && x < 50) {if(y>49) y=49; correction =  _5B70[y] - _5B70[x] ;lbe=  _5B70[y];correctL=true;}    //

					else if( (hue >= (_15PB70[x] - 0.035)) && (hue < (_15PB70[x] + _3PB70[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _15PB70[y] - _15PB70[x] ;lbe=_15PB70[y];correctL=true; }
					else if (( hue>=(_45PB70[x] + _3PB70[x])/2.0)  && (hue < (_6PB70[x] + _45PB70[x])/2.0) && x < 50) {if(y>49) y=49;correction =  _45PB70[y] - _45PB70[x] ;lbe=_45PB70[y];correctL=true;}
					else if (( hue>=(_6PB70[x] + _45PB70[x])/2.0)  && (hue < (_6PB70[x] + _75PB70[x])/2.0) && x <50) {if(y>49) y=49;correction =  _6PB70[y] - _6PB70[x] ;lbe=_6PB70[y];correctL=true;}
					else if (( hue>=(_6PB70[x] + _75PB70[x])/2.0)  && (hue < (_9PB70[x] + _75PB70[x])/2.0) && x <50) {if(y>49) y=49;correction =  _75PB70[y] - _75PB70[x] ;lbe=_75PB70[y];correctL=true; }
					else if (( hue>=(_9PB70[x] + _75PB70[x])/2.0)  && (hue < (_9PB70[x] + 0.035)) && x <50) {if(y>49) y=49;correction =  _9PB70[y] - _9PB70[x] ;lbe=_9PB70[y];correctL=true;}
				}
				else if (lum < 85.0) {
					if( (hue <= (_05PB80[x] + _15PB80[x])/2.0) && (hue > (_05PB80[x] + _10B80[x])/2.0) && x<40 ) {if(y>39) y=39; correction =  _05PB80[y] - _05PB80[x] ;lbe=_05PB80[y] ;correctL=true;}
					else if( (hue <= (_05PB80[x] + _10B80[x])/2.0) && (hue >(_10B80[x] + _9B80[x])/2.0) && x<40 ) {if(y>39) y=39;correction =  _10B80[y] - _10B80[x] ;lbe=_10B80[y];correctL=true;}
					else if( (hue <= (_10B80[x] + _9B80[x])/2.0) && (hue >(_9B80[x] + _7B80[x])/2.0) && x<40 ) {if(y>39) y=39;correction =  _9B80[y] - _9B80[x] ;lbe= _9B80[y];correctL=true;}
					else if( (hue <= (_9B80[x] + _7B80[x])/2.0) && (hue >(_5B80[x] + _7B80[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _7B80[y] - _7B80[x] ;lbe=_7B80[y];correctL=true;}
					else if (( hue<=(_5B80[x] + _7B80[x])/2.0)  && (hue > (_5B80[x]-0.035)) && x < 50) {if(y>49) y=49; correction =  _5B80[y] - _5B80[x] ; lbe=_5B80[y];correctL=true;}    //

					else if( (hue >= (_15PB80[x] - 0.035)) && (hue < (_15PB80[x] + _3PB80[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _15PB80[y] - _15PB80[x] ; lbe=_15PB80[y];correctL=true;}
					else if (( hue>=(_45PB80[x] + _3PB80[x])/2.0)  && (hue < (_6PB80[x] + _45PB80[x])/2.0) && x < 50) {if(y>49) y=49;correction =  _45PB80[y] - _45PB80[x] ;lbe= _45PB80[y];correctL=true;}
					else if (( hue>=(_6PB80[x] + _45PB80[x])/2.0)  && (hue < (_6PB80[x] + _75PB80[x])/2.0) && x <50) {if(y>49) y=49;correction =  _6PB80[y] - _6PB80[x] ;lbe=_6PB80[y];correctL=true;}
					else if (( hue>=(_6PB80[x] + _75PB80[x])/2.0)  && (hue < (_9PB80[x] + _75PB80[x])/2.0) && x <50) {if(y>49) y=49;correction =  _75PB80[y] - _75PB80[x] ;lbe=_75PB80[y];correctL=true; }
					else if (( hue>=(_9PB80[x] + _75PB80[x])/2.0)  && (hue < (_9PB80[x] + 0.035)) && x <50) {if(y>49) y=49;correction =  _9PB80[y] - _9PB80[x] ;lbe=_9PB80[y]; correctL=true;}
				}
            }
        }
        // end PB correction

        //red yellow correction
        else if(zone==2) {
            if(lum > 15.0) {
				if(lum < 25.0) {
					if( (hue <= (_10YR20[x] + 0.035)) && (hue > (_10YR20[x] + _85YR20[x])/2.0) && x<=45) {if(y>49) y=49;correction =  _10YR20[y] - _10YR20[x] ;lbe=_10YR20[y];correctL=true;}
					else if (( hue<=(_85YR20[x] + _10YR20[x])/2.0)  && (hue > (_85YR20[x] + 0.035) && x <= 45)) {if(y>49) y=49;correction =  _85YR20[y] - _85YR20[x] ;lbe= _85YR20[y];correctL=true;}
				}
				else if (lum <35.0) {
					if( (hue <= (_10YR30[x] + 0.035)) && (hue > (_10YR30[x] + _85YR30[x])/2.0) && x < 85) {if(y>89) y=89;correction =  _10YR30[y] - _10YR30[x] ;lbe=_10YR30[y];correctL=true;}
					else if( (hue <= (_10YR30[x] + _85YR30[x])/2.0) && (hue >(_85YR30[x] + _7YR30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _85YR30[y] - _85YR30[x] ;lbe= _85YR30[y];correctL=true;}
					else if (( hue<=(_85YR30[x] + _7YR30[x])/2.0)  && (hue > (_7YR30[x] + _55YR30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _7YR30[y] - _7YR30[x] ;lbe=_7YR30[y];correctL=true;}
					else if (( hue<=(_7YR30[x] + _55YR30[x])/2.0)  && (hue > (_55YR30[x] + _4YR30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _55YR30[y] - _55YR30[x] ;lbe=_55YR30[y];correctL=true; }
					else if (( hue<=(_55YR30[x] + _4YR30[x])/2.0)  && (hue > (_4YR30[x] + _25YR30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _4YR30[y] - _4YR30[x] ;lbe= _4YR30[y];correctL=true;}
					else if (( hue<=(_4YR30[x] + _25YR30[x])/2.0)  && (hue > (_25YR30[x] + _10R30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR30[y] - _25YR30[x] ;lbe=_25YR30[y];correctL=true;}
					else if (( hue<=(_25YR30[x] + _10R30[x])/2.0)  && (hue > (_10R30[x] + _9R30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10R30[y] - _10R30[x] ; lbe=_10R30[y];correctL=true;}
					else if (( hue<=(_10R30[x] + _9R30[x])/2.0)  && (hue > (_9R30[x] + _7R30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _9R30[y] - _9R30[x] ;lbe=_9R30[y];correctL=true;}
					else if (( hue<=(_9R30[x] + _7R30[x])/2.0)  && (hue > (_7R30[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _7R30[y] - _7R30[x] ; lbe=_7R30[y] ;correctL=true;}
				}
				else if (lum <45.0) {
					if( (hue <= (_10YR40[x] + 0.035)) && (hue > (_10YR40[x] + _85YR40[x])/2.0)&& x<85) {if(y>89) y=89;correction =  _10YR40[y] - _10YR40[x] ;lbe=_10YR40[y];correctL=true;}
					else if( (hue <= (_10YR40[x] + _85YR40[x])/2.0) && (hue >(_85YR40[x] + _7YR40[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _85YR40[y] - _85YR40[x] ;lbe= _85YR40[y];correctL=true;}
					else if (( hue<=(_85YR40[x] + _7YR40[x])/2.0)  && (hue > (_7YR40[x] + _55YR40[x])/2.0) && x < 85) {if(y>89) y=89;correction =  _7YR40[y] - _7YR40[x] ;lbe= _7YR40[y];correctL=true;}
					else if (( hue<=(_7YR40[x] + _55YR40[x])/2.0)  && (hue > (_55YR40[x] + _4YR40[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _55YR40[y] - _55YR40[x] ;lbe=_55YR40[y];correctL=true; }
					else if (( hue<=(_55YR40[x] + _4YR40[x])/2.0)  && (hue > (_4YR40[x] + _25YR40[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _4YR40[y] - _4YR40[x] ;lbe=_4YR40[y]; correctL=true;}
					else if (( hue<=(_4YR40[x] + _25YR40[x])/2.0)  && (hue > (_25YR40[x] + _10R40[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR40[y] - _25YR40[x] ;lbe=_25YR40[y] ;correctL=true;}
					else if (( hue<=(_25YR40[x] + _10R40[x])/2.0)  && (hue > (_10R40[x] + _9R40[x])/2.0) && x < 85) {if(y>89) y=89;correction =  _10R40[y] - _10R40[x] ; lbe=_10R40[y];correctL=true;}
					else if (( hue<=(_10R40[x] + _9R40[x])/2.0)  && (hue > (_9R40[x] + _7R40[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _9R40[y] - _9R40[x] ;lbe=_9R40[y];correctL=true;}
					else if (( hue<=(_9R40[x] + _7R40[x])/2.0)  && (hue > (_7R40[x] -0.035))&& x < 85 ) {if(y>89) y=89;correction =  _7R40[y] - _7R40[x] ; lbe=_7R40[y];correctL=true;}
				}
				else if (lum <55.0) {
					if( (hue <= (_10YR50[x] + 0.035)) && (hue > (_10YR50[x] + _85YR50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10YR50[y] - _10YR50[x] ;lbe=_10YR50[y];correctL=true;}
					else if( (hue <= (_10YR50[x] + _85YR50[x])/2.0) && (hue >(_85YR50[x] + _7YR50[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _85YR50[y] - _85YR50[x] ;lbe=_85YR50[y];correctL=true;}
					else if (( hue<=(_85YR50[x] + _7YR50[x])/2.0)  && (hue > (_7YR50[x] + _55YR50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _7YR50[y] - _7YR50[x] ;lbe=_7YR50[y];correctL=true;}
					else if (( hue<=(_7YR50[x] + _55YR50[x])/2.0)  && (hue > (_55YR50[x] + _4YR50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _55YR50[y] - _55YR50[x] ; lbe=_55YR50[y];correctL=true;}
					else if (( hue<=(_55YR50[x] + _4YR50[x])/2.0)  && (hue > (_4YR50[x] + _25YR50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _4YR50[y] - _4YR50[x] ;lbe=_4YR50[y]; correctL=true;}
					else if (( hue<=(_4YR50[x] + _25YR50[x])/2.0)  && (hue > (_25YR50[x] + _10R50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR50[y] - _25YR50[x] ;lbe=_25YR50[y];correctL=true;}
					else if (( hue<=(_25YR50[x] + _10R50[x])/2.0)  && (hue > (_10R50[x] + _9R50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10R50[y] - _10R50[x] ;lbe= _10R50[y]; correctL=true;}
					else if (( hue<=(_10R50[x] + _9R50[x])/2.0)  && (hue > (_9R50[x] + _7R50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _9R50[y] - _9R50[x] ;lbe=_9R50[y];correctL=true;}
					else if (( hue<=(_9R50[x] + _7R50[x])/2.0)  && (hue > (_7R50[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _7R50[y] - _7R50[x] ; lbe=_7R50[y];correctL=true;}
				}
				else if (lum <65.0) {
					if( (hue <= (_10YR60[x] + 0.035)) && (hue > (_10YR60[x] + _85YR60[x])/2.0)) {;correction =  _10YR60[y] - _10YR60[x] ;lbe= _10YR60[y];correctL=true;}
					else if( (hue <= (_10YR60[x] + _85YR60[x])/2.0) && (hue >(_85YR60[x] + _7YR60[x])/2.0) ) {;correction =  _85YR60[y] - _85YR60[x] ;lbe= _85YR60[y];correctL=true;}
					else if (( hue<=(_85YR60[x] + _7YR60[x])/2.0)  && (hue > (_7YR60[x] + _55YR60[x])/2.0)) {correction =  _7YR60[y] - _7YR60[x] ;lbe=_7YR60[y];correctL=true;}
					else if (( hue<=(_7YR60[x] + _55YR60[x])/2.0)  && (hue > (_55YR60[x] + _4YR60[x])/2.0)) {correction =  _55YR60[y] - _55YR60[x] ;lbe= _55YR60[y];correctL=true;}
					else if (( hue<=(_55YR60[x] + _4YR60[x])/2.0)  && (hue > (_4YR60[x] + _25YR60[x])/2.0)) {correction =  _4YR60[y] - _4YR60[x] ;lbe=_4YR60[y]; correctL=true;}
					else if (( hue<=(_4YR60[x] + _25YR60[x])/2.0)  && (hue > (_25YR60[x] + _10R60[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR60[y] - _25YR60[x] ;lbe=_25YR60[y];correctL=true;}
					else if (( hue<=(_25YR60[x] + _10R60[x])/2.0)  && (hue > (_10R60[x] + _9R60[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10R60[y] - _10R60[x] ;lbe= _10R60[y];correctL=true; }
					else if (( hue<=(_10R60[x] + _9R60[x])/2.0)  && (hue > (_9R60[x] + _7R60[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _9R60[y] - _9R60[x] ;lbe=_9R60[y];correctL=true;}
					else if (( hue<=(_9R60[x] + _7R60[x])/2.0)  && (hue > (_7R60[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _7R60[y] - _7R60[x] ;lbe=_7R60[y]; correctL=true;}
				}
				else if (lum <75.0) {
					if( (hue <= (_10YR70[x] + 0.035)) && (hue > (_10YR70[x] + _85YR70[x])/2.0)) {correction =  _10YR70[y] - _10YR70[x] ;lbe= _10YR70[y];correctL=true;}
					else if( (hue <= (_10YR70[x] + _85YR70[x])/2.0) && (hue >(_85YR70[x] + _7YR70[x])/2.0)) {correction =  _85YR70[y] - _85YR70[x] ;lbe=_85YR70[y];correctL=true;}
					 if (( hue<=(_85YR70[x] + _7YR70[x])/2.0)  && (hue > (_7YR70[x] + _55YR70[x])/2.0)) {correction =  _7YR70[y] - _7YR70[x] ;lbe=_7YR70[y];correctL=true;}
					else if (( hue<=(_7YR70[x] + _55YR70[x])/2.0)  && (hue > (_55YR70[x] + _4YR70[x])/2.0)) {correction =  _55YR70[y] - _55YR70[x] ;lbe=_55YR70[y];correctL=true; }
					else if (( hue<=(_55YR70[x] + _4YR70[x])/2.0)  && (hue > (_4YR70[x] + _25YR70[x])/2.0)) {correction =  _4YR70[y] - _4YR70[x] ;lbe=_4YR70[y];correctL=true; }
					else if (( hue<=(_4YR70[x] + _25YR70[x])/2.0)  && (hue > (_25YR70[x] + _10R70[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR70[y] - _25YR70[x] ;lbe= _25YR70[y];correctL=true;}
					else if (( hue<=(_25YR70[x] + _10R70[x])/2.0)  && (hue > (_10R70[x] + _9R70[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10R70[y] - _10R70[x] ;lbe= _10R70[y];correctL=true;}
					else if (( hue<=(_10R70[x] + _9R70[x])/2.0)  && (hue > (_9R70[x] + _7R70[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _9R70[y] - _9R70[x] ;lbe= _9R70[y] ;correctL=true;}
					else if (( hue<=(_9R70[x] + _7R70[x])/2.0)  && (hue > (_7R70[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _7R70[y] - _7R70[x] ;lbe=_7R70[y];correctL=true; }
				}
				else if (lum <85.0) {
					if( (hue <= (_10YR80[x] + 0.035)) && (hue > (_10YR80[x] + _85YR80[x])/2.0)) {correction =  _10YR80[y] - _10YR80[x] ;lbe=_10YR80[y];correctL=true;}
					else if( (hue <= (_10YR80[x] + _85YR80[x])/2.0) && (hue >(_85YR80[x] + _7YR80[x])/2.0)) {correction =  _85YR80[y] - _85YR80[x] ;lbe= _85YR80[y];}
					else if (( hue<=(_85YR80[x] + _7YR80[x])/2.0)  && (hue > (_7YR80[x] + _55YR80[x])/2.0) && x<85) {if(y>89) y=89;correction =  _7YR80[y] - _7YR80[x] ;lbe=_7YR80[y];correctL=true;}
					else if (( hue<=(_7YR80[x] + _55YR80[x])/2.0)  && (hue > (_55YR80[x] + _4YR80[x])/2.0) && x <45) {correction =  _55YR80[y] - _55YR80[x] ;lbe=_55YR80[y];correctL=true; }
					else if (( hue<=(_55YR80[x] + _4YR80[x])/2.0)  && (hue > (_4YR80[x] - 0.035) && x<45)) {if(y>49) y=49;correction =  _4YR80[y] - _4YR80[x] ; lbe=_4YR80[y] ;correctL=true;}
				}
				else if (lum <95.0) {
					if( (hue <= (_10YR90[x] + 0.035)) && (hue > (_10YR90[x] -0.035) && x<85)) {if(y>89) y=89;correction =  _10YR90[y] - _10YR90[x] ;lbe= _10YR90[y];correctL=true;}
					else if ( hue<=(_85YR90[x] + 0.035)  && hue > (_85YR90[x] -0.035) && x<85) {if(y>89) y=89;correction =  _85YR90[y] - _85YR90[x] ;lbe=_85YR90[y];correctL=true;}
					else if (( hue<=(_55YR90[x] + 0.035)  && (hue > (_55YR90[x] - 0.035) && x<45))) {if(y>49) y=49;correction =  _55YR90[y] - _55YR90[x] ;lbe= _55YR90[y];correctL=true; }
				}
            }
        }
        //end red yellow

        //Green yellow correction
        else if(zone==3) {
            if (lum >= 25.0) {
				if (lum <35.0) {
					if( (hue <= (_7G30[x] + 0.035)) && (hue > (_7G30[x] + _5G30[x])/2.0) ) {correction =  _7G30[y] - _7G30[x] ;lbe=_7G30[y];correctL=true;}
					else if( (hue <= (_7G30[x] + _5G30[x])/2.0) && (hue >(_5G30[x] + _25G30[x])/2.0)) {correction =  _5G30[y] - _5G30[x] ;lbe= _5G30[y];correctL=true;}
					else if (( hue<=(_25G30[x] + _5G30[x])/2.0)  && (hue > (_25G30[x] + _1G30[x])/2.0)) {correction =  _25G30[y] - _25G30[x] ;lbe=_25G30[y];correctL=true;}
					else if (( hue<=(_1G30[x] + _25G30[x])/2.0)  && (hue > (_1G30[x] + _10GY30[x])/2.0)) {correction =  _1G30[y] - _1G30[x] ;lbe= _1G30[y];correctL=true;}
					else if (( hue<=(_1G30[x] + _10GY30[x])/2.0)  && (hue > (_10GY30[x] + _75GY30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10GY30[y] - _10GY30[x] ;lbe=  _10GY30[y];correctL=true;}
					else if (( hue<=(_10GY30[x] + _75GY30[x])/2.0)  && (hue > (_75GY30[x] + _5GY30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _75GY30[y] - _75GY30[x] ;lbe=_75GY30[y];correctL=true;}
					else if (( hue<=(_5GY30[x] + _75GY30[x])/2.0)  && (hue > (_5GY30[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _5GY30[y] - _5GY30[x] ;lbe= _5GY30[y] ;correctL=true; }
				}
				else if (lum <45.0) {
					if( (hue <= (_7G40[x] + 0.035)) && (hue > (_7G40[x] + _5G40[x])/2.0) ) {correction =  _7G40[y] - _7G40[x] ;lbe= _7G40[y];correctL=true;}
					else if( (hue <= (_7G40[x] + _5G40[x])/2.0) && (hue >(_5G40[x] + _25G40[x])/2.0)) {correction =  _5G40[y] - _5G40[x] ;lbe=_5G40[y];correctL=true;}
					else if (( hue<=(_25G40[x] + _5G40[x])/2.0)  && (hue > (_25G40[x] + _1G40[x])/2.0)) {correction =  _25G40[y] - _25G40[x] ;lbe=_25G40[y];correctL=true;}
					else if (( hue<=(_1G40[x] + _25G40[x])/2.0)  && (hue > (_1G40[x] + _10GY40[x])/2.0)) {correction =  _1G40[y] - _1G40[x] ;lbe=_1G40[y];correctL=true; }
					else if (( hue<=(_1G40[x] + _10GY40[x])/2.0)  && (hue > (_10GY40[x] + _75GY40[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10GY40[y] - _10GY40[x] ;lbe=_10GY40[y];correctL=true; }
					else if (( hue<=(_10GY40[x] + _75GY40[x])/2.0)  && (hue > (_75GY40[x] + _5GY40[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _75GY40[y] - _75GY40[x] ;lbe=_75GY40[y];correctL=true;}
					else if (( hue<=(_5GY40[x] + _75GY40[x])/2.0)  && (hue > (_5GY40[x]-0.035)) && x < 85) {if(y>89) y=89; correction =  _5GY40[y] - _5GY40[x] ;lbe=_5GY40[y];correctL=true; }    //
				}
				else if (lum <55.0) {
					if( (hue <= (_7G50[x] + 0.035)) && (hue > (_7G50[x] + _5G50[x])/2.0) ) {correction =  _7G50[y] - _7G50[x] ;lbe=_7G50[y];correctL=true;}
					else if( (hue <= (_7G50[x] + _5G50[x])/2.0) && (hue >(_5G50[x] + _25G50[x])/2.0)) {correction =  _5G50[y] - _5G50[x] ;lbe=_5G50[y];correctL=true;}
					else if (( hue<=(_25G50[x] + _5G50[x])/2.0)  && (hue > (_25G50[x] + _1G50[x])/2.0)) {correction =  _25G50[y] - _25G50[x] ;lbe= _25G50[y];correctL=true;}
					else if (( hue<=(_1G50[x] + _25G50[x])/2.0)  && (hue > (_1G50[x] + _10GY50[x])/2.0)) {correction =  _1G50[y] - _1G50[x] ; lbe=_1G50[y];correctL=true;}
					else if (( hue<=(_1G50[x] + _10GY50[x])/2.0)  && (hue > (_10GY50[x] + _75GY50[x])/2.0)) {correction =  _10GY50[y] - _10GY50[x] ;lbe= _10GY50[y];correctL=true;}
					else if (( hue<=(_10GY50[x] + _75GY50[x])/2.0)  && (hue > (_75GY50[x] + _5GY50[x])/2.0)) {correction =  _75GY50[y] - _75GY50[x] ;lbe=_75GY50[y];correctL=true;}
					else if (( hue<=(_5GY50[x] + _75GY50[x])/2.0)  && (hue > (_5GY50[x] -0.035))) {correction =  _5GY50[y] - _5GY50[x] ; lbe=_5GY50[y];correctL=true;}
				}
				else if (lum <65.0) {
					if( (hue <= (_7G60[x] + 0.035)) && (hue > (_7G60[x] + _5G60[x])/2.0) ) {correction =  _7G60[y] - _7G60[x] ;lbe=_7G60[y];correctL=true;}
					else if( (hue <= (_7G60[x] + _5G60[x])/2.0) && (hue >(_5G60[x] + _25G60[x])/2.0)) {correction =  _5G60[y] - _5G60[x] ;lbe=_5G60[y];correctL=true;}
					else if (( hue<=(_25G60[x] + _5G60[x])/2.0)  && (hue > (_25G60[x] + _1G60[x])/2.0)) {correction =  _25G60[y] - _25G60[x] ;lbe=_25G60[y];correctL=true;}
					else if (( hue<=(_1G60[x] + _25G60[x])/2.0)  && (hue > (_1G60[x] + _10GY60[x])/2.0)) {correction =  _1G60[y] - _1G60[x] ; lbe=_1G60[y];correctL=true;}
					else if (( hue<=(_1G60[x] + _10GY60[x])/2.0)  && (hue > (_10GY60[x] + _75GY60[x])/2.0)) {correction =  _10GY60[y] - _10GY60[x] ;lbe= _10GY60[y];correctL=true;}
					else if (( hue<=(_10GY60[x] + _75GY60[x])/2.0)  && (hue > (_75GY60[x] + _5GY60[x])/2.0)) {correction =  _75GY60[y] - _75GY60[x] ;lbe=_75GY60[y] ;correctL=true;}
					else if (( hue<=(_5GY60[x] + _75GY60[x])/2.0)  && (hue > (_5GY60[x] -0.035))) {correction =  _5GY60[y] - _5GY60[x] ; lbe=_5GY60[y];correctL=true;}
				}
				else if (lum <75.0) {
					if( (hue <= (_7G70[x] + 0.035)) && (hue > (_7G70[x] + _5G70[x])/2.0) ) {correction =  _7G70[y] - _7G70[x] ;lbe= _7G70[y];correctL=true;}
					else if( (hue <= (_7G70[x] + _5G70[x])/2.0) && (hue >(_5G70[x] + _25G70[x])/2.0)) {correction =  _5G70[y] - _5G70[x] ;lbe= _5G70[y];correctL=true;}
					else if (( hue<=(_25G70[x] + _5G70[x])/2.0)  && (hue > (_25G70[x] + _1G70[x])/2.0)) {correction =  _25G70[y] - _25G70[x] ;lbe=_25G70[y];correctL=true;}
					else if (( hue<=(_1G70[x] + _25G70[x])/2.0)  && (hue > (_1G70[x] + _10GY70[x])/2.0)) {correction =  _1G70[y] - _1G70[x] ;lbe= _1G70[y] ;correctL=true;}
					else if (( hue<=(_1G70[x] + _10GY70[x])/2.0)  && (hue > (_10GY70[x] + _75GY70[x])/2.0)) {correction =  _10GY70[y] - _10GY70[x] ; lbe=_10GY70[y];correctL=true;}
					else if (( hue<=(_10GY70[x] + _75GY70[x])/2.0)  && (hue > (_75GY70[x] + _5GY70[x])/2.0)) {correction =  _75GY70[y] - _75GY70[x] ;lbe=_75GY70[y];correctL=true;}
					else if (( hue<=(_5GY70[x] + _75GY70[x])/2.0)  && (hue > (_5GY70[x] -0.035))) {correction =  _5GY70[y] - _5GY70[x] ;lbe=  _5GY70[y];correctL=true;}
				}
				else if (lum <85.0) {
					if( (hue <= (_7G80[x] + 0.035)) && (hue > (_7G80[x] + _5G80[x])/2.0) ) {correction =  _7G80[y] - _7G80[x] ;lbe=_7G80[y];correctL=true;}
					else if( (hue <= (_7G80[x] + _5G80[x])/2.0) && (hue >(_5G80[x] + _25G80[x])/2.0)) {correction =  _5G80[y] - _5G80[x] ;lbe=_5G80[y];correctL=true;}
					else if (( hue<=(_25G80[x] + _5G80[x])/2.0)  && (hue > (_25G80[x] + _1G80[x])/2.0)) {correction =  _25G80[y] - _25G80[x] ;lbe=_25G80[y];correctL=true;}
					else if (( hue<=(_1G80[x] + _25G80[x])/2.0)  && (hue > (_1G80[x] + _10GY80[x])/2.0)) {correction =  _1G80[y] - _1G80[x] ; lbe= _1G80[y];correctL=true;}
					else if (( hue<=(_1G80[x] + _10GY80[x])/2.0)  && (hue > (_10GY80[x] + _75GY80[x])/2.0)) {correction =  _10GY80[y] - _10GY80[x] ;lbe=_10GY80[y];correctL=true; }
					else if (( hue<=(_10GY80[x] + _75GY80[x])/2.0)  && (hue > (_75GY80[x] + _5GY80[x])/2.0)) {correction =  _75GY80[y] - _75GY80[x] ;lbe=_75GY80[y];correctL=true;}
					else if (( hue<=(_5GY80[x] + _75GY80[x])/2.0)  && (hue > (_5GY80[x] -0.035))) {correction =  _5GY80[y] - _5GY80[x] ; lbe=_5GY80[y];correctL=true;}
				}
            }
        }
        //end green yellow

        //Red purple correction : only for L < 30
        else if(zone==4) {
            if (lum > 5.0) {
				if (lum < 15.0) {
					if( (hue <= (_5R10[x] + 0.035)) && (hue > (_5R10[x] - 0.043)) && x<45) {if(y>44) y=44;correction =  _5R10[y] - _5R10[x] ;lbe=_5R10[y];correctL=true;}
					else if( (hue <= (_25R10[x] + 0.043)) && (hue >(_25R10[x] + _10RP10[x])/2.0) && x<45 ) {if(y>44) y=44;correction =  _25R10[y] - _25R10[x] ;lbe= _25R10[y];correctL=true;}
					else if ( (hue <=(_25R10[x] + _10RP10[x])/2.0) && (hue > (_10RP10[x] -0.035) ) && x<45){if(y>44) y=44; correction =  _10RP10[y] - _10RP10[x] ;lbe=_10RP10[y];correctL=true;}
				}
				else if (lum <25.0) {
					if( (hue <= (_5R20[x] + 0.035)) && (hue > (_5R20[x] + _25R20[x])/2.0) && x<70 ) {if(y>70) y=70;correction =  _5R20[y] - _5R20[x] ;lbe= _5R20[y];correctL=true;}
					else if( (hue <= (_5R20[x] + _25R20[x])/2.0) && (hue >(_10RP20[x] + _25R20[x])/2.0) && x<70) {if(y>70) y=70;correction =  _25R20[y] - _25R20[x] ;lbe=_25R20[y];correctL=true;}
					else if (( hue<=(_10RP20[x] + _25R20[x])/2.0)  && (hue > (_10RP20[x] -0.035)) && x<70) {if(y>70) y=70; correction =  _10RP20[y] - _10RP20[x] ;lbe=_10RP20[y];correctL=true;}
				}
				else if (lum <35.0) {
					if( (hue <= (_5R30[x] + 0.035)) && (hue > (_5R30[x] + _25R30[x])/2.0) && x<85 ) {if(y>85) y=85;correction =  _5R30[y] - _5R30[x] ;lbe= _5R30[y];correctL=true;}
					else if( (hue <= (_5R30[x] + _25R30[x])/2.0) && (hue >(_10RP30[x] + _25R30[x])/2.0) && x< 85) {if(y>85) y=85;correction =  _25R30[y] - _25R30[x] ;lbe=_25R30[y];correctL=true;}
					else if (( hue<=(_10RP30[x] + _25R30[x])/2.0)  && (hue > (_10RP30[x] -0.035)) && x<85) {if(y>85) y=85; correction =  _10RP30[y] - _10RP30[x] ;lbe= _10RP30[y];correctL=true;}
				}
            }
        }
        //end red purple
    }


    /*
     * SkinSat
     * Copyright (c)2011  Jacques Desmis <jdesmis@gmail.com>
     *
     * skin color: mixed from NX2 skin color palette, Von Luschan, and photos of people white,
     * black, yellow....there are some little exceptions...cover 99% case
     * pay attention to white balance, and do not change hue and saturation, upstream of the modification
     *
     */
    void Color::SkinSat (float lum, float hue, float chrom, float &satreduc, int chromx) {

        float reduction=0.3;// to be adapted...by tests
        float extendedreduction=0.4;
        float extendedreduction2=0.6;

        if(chromx==1) {reduction=0.6;extendedreduction=0.7;extendedreduction2=0.8;}

        float C9=0.0, C8=0.0, C7=0.0, C4=0.0, C3=0.0, C2=0.0, C1=0.0;
        float H9=0.0, H8=0.0, H7=0.0, H4=0.0, H3=0.0, H2=0.0, H1=0.0, H10=0.0,H11=0.0;
        H9=0.05;H8=0.25;H7=0.1;H4=0.02;H3=0.02;H2=0.1;H1=0.1;H10=-0.2;H11=-0.2;//H10 and H11 are curious...H11=-0.8 ??
        C9=8.0;C8=15.0;C7=12.0;C4=7.0;C3=5.0;C2=5.0;C1=5.0;

        // wide area for transition

        if      (lum >= 92.0 && (hue > -0.1 && hue < 1.65) && (chrom > 7.0 && chrom < (18.0))) satreduc=extendedreduction2;
        else if (lum >= 85.0 && lum < 92.0 && (hue > 0.0 && hue < 1.65) && (chrom > 7.0 && chrom < (35.0+C9))) satreduc=extendedreduction2;
        else if ((lum > 20 && lum < 85) && (hue > (0.02 + H11) && hue < 1.65) && (chrom > 7.0 && chrom < (55.0+C9) )) satreduc=extendedreduction2;
        else if (lum < 20.0 && (hue > (0.02+H11) && hue < 1.60) && (chrom > 7.0 && chrom < (45.0+C1) )) satreduc=extendedreduction2;

        // wide area  skin color, useful if not accurate colorimetry or if the user has changed hue and saturation

        if      (lum >= 92.0  && (hue > 0.8 && hue < 1.65) && (chrom > 7.0 && chrom < (15.0))) satreduc=extendedreduction;
        else if (lum >= 85.0 && lum < 92.0  && (hue > 0.70 && hue < 1.4) && (chrom > 7.0 && chrom < (26.0+C9))) satreduc=extendedreduction;
        else if ((lum > 20 && lum < 85) && (hue > (0.02 + H11) && hue < 1.5) && (chrom > 7.0 && chrom < (48.0+C9) )) satreduc=extendedreduction;
        else if (lum < 20.0  && (hue > (0.02+H11) && hue < 1.0) && (chrom > 7.0 && chrom < (35.0+C1) )) satreduc=extendedreduction;

        // "real" skin color : take into account a slightly usage of contrast and saturation in RT if option "skin" = 1

        if       (lum >= 85.0  && (hue > (0.78-H9) && hue < (1.18+H9)) && (chrom > 8.0 && chrom < (14.0+C9))) satreduc=reduction;
        else if ((lum >= 70.0 && lum < 85.0)  && (hue > 0.4 && hue < (1.04+H8)) && (chrom > 8.0 && chrom < (35.0+C8))) satreduc=reduction;
        else if ((lum >= 52.0 && lum < 70.0)  && (hue > 0.3 && hue < (1.27+H7)) && (chrom > 11.0 && chrom < (35.0+C7))) satreduc=reduction;
        else if ((lum >= 35.0 && lum < 52.0)  && (hue > 0.3 && hue < (1.25+H4)) && (chrom > 13.0 && chrom < (37.0+C4))) satreduc=reduction;
        else if ((lum >= 20.0 && lum < 35.0)  && (hue > 0.3 && hue < (1.20+H3)) && (chrom > 7.0 && chrom <(35.0+C3) )) satreduc=reduction;
        else if ((lum > 10.0 && lum < 20.0)  && (hue > (0.0 + H10) && hue < (0.95 +H2)) && (chrom > 8.0 && chrom < (23.0+C2))) satreduc=reduction;
        else if ((lum < 10.0)  && (hue > (0.02 + H10) && hue < (0.90+H1)) && (chrom > 8.0 && chrom < (23.0+C1))) satreduc=reduction; // no data : extrapolate
    }

    /*
     * Munsell Lch correction
     * Copyright (c) 2011  Jacques Desmis <jdesmis@gmail.com>
     *
     * data (Munsell ==> Lab) obtained with WallKillcolor and http://www.cis.rit.edu/research/mcsl2/online/munsell.php
     * each LUT give Hue in function of C, for each color Munsell and Luminance
     * eg: _6PB20 : color Munsell 6PB for L=20 c=5 c=45 c=85 c=125..139 when possible: interpolation betwwen values
     * no value for C<5  (gray)
     * low memory footprint -- maximum: 195 LUTf * 140 values
     * errors due to small number of samples in LUT and linearization are very low (1 to 2%)
     * errors due to a different illuminant "Daylight" than "C" are low, about 10%. For example, a theoretical correction of 0.1 radian will be made with a real correction of 0.09 or 0.11 depending on the color illuminant D50
     * errors due to the use of a very different illuminant "C", for example illuminant "A" (tungsten) are higher, about 20%. Theoretical correction of 0.52 radians will be made with a real correction of 0.42
     */
    void Color::initMunsell () {
    #ifdef _DEBUG
        MyTime t1e,t2e;
        t1e.set();
    #endif

    const int maxInd  = 140;
    const int maxInd2 = 90;
    const int maxInd3 = 50;

    //blue for sky
    _5B40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _5B40[i] = -2.3 + 0.0025*(i-5);
            else if (i<90 && i>=45) _5B40[i] = -2.2 + 0.00*(i-45);
        }
        //printf("5B %1.2f  %1.2f\n",_5B40[44],_5B40[89]);
    _5B50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _5B50[i] = -2.34 + 0.0025*(i-5);
            else if (i<90 && i>=45) _5B50[i] = -2.24+0.0003*(i-45);
        }
        //printf("5B %1.2f  %1.2f\n",_5B50[44],_5B50[89]);
    _5B60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _5B60[i] = -2.4 + 0.003*(i-5);
            else if (i<90 && i>=45) _5B60[i] = -2.28+0.0005*(i-45);
        }
        //printf("5B %1.2f  %1.2f\n",_5B60[44],_5B60[89]);
    _5B70(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _5B70[i] = -2.41 + 0.00275*(i-5);
            else if (i<90 && i>=45) _5B70[i] = -2.30+0.00025*(i-45);
        }
        //printf("5B %1.2f  %1.2f\n",_5B70[44],_5B70[89]);
    _5B80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _5B80[i] = -2.45 +0.003*(i-5);
        }
        //printf("5B %1.2f\n",_5B80[49]);

    _7B40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7B40[i] = -2.15 + 0.0027*(i-5);
            else if (i<90 && i>=45) _7B40[i] = -2.04 + 0.00*(i-45);
        }
        //printf("7B %1.2f  %1.2f\n",_7B40[44],_7B40[89]);
    _7B50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7B50[i] = -2.20 + 0.003*(i-5);
            else if (i<90 && i>=45) _7B50[i] = -2.08 + 0.001*(i-45);
        }
        //printf("7B %1.2f  %1.2f\n",_7B50[44],_7B50[79]);
    _7B60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7B60[i] = -2.26 + 0.0035*(i-5);
            else if (i<90 && i>=45) _7B60[i] = -2.12 + 0.001*(i-45);
        }
        //printf("7B %1.2f  %1.2f\n",_7B60[44],_7B60[79]);
    _7B70(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7B70[i] = -2.28 + 0.003*(i-5);
            else if (i<90 && i>=45) _7B70[i] = -2.16 + 0.0015*(i-45);
        }
        //printf("7B %1.2f  %1.2f\n",_7B70[44],_7B70[64]);
    _7B80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _7B80[i] = -2.30 +0.0028*(i-5);
        }
        //printf("5B %1.2f\n",_7B80[49]);

    _9B40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _9B40[i] = -1.99 + 0.0022*(i-5);
            else if (i<90 && i>=45) _9B40[i] = -1.90 + 0.0008*(i-45);
        }
        //printf("9B %1.2f  %1.2f\n",_9B40[44],_9B40[69]);
    _9B50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _9B50[i] = -2.04 + 0.0025*(i-5);
            else if (i<90 && i>=45) _9B50[i] = -1.94 + 0.0013*(i-45);
        }
        //printf("9B %1.2f  %1.2f\n",_9B50[44],_9B50[77]);
    _9B60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _9B60[i] = -2.10 + 0.0033*(i-5);
            else if (i<90 && i>=45) _9B60[i] = -1.97 + 0.001*(i-45);
        }
        //printf("9B %1.2f  %1.2f\n",_9B60[44],_9B60[79]);
    _9B70(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _9B70[i] = -2.12 + 0.003*(i-5);
            else if (i<90 && i>=45) _9B70[i] = -2.00 + 0.001*(i-45);
        }
        //printf("9B %1.2f  %1.2f\n",_9B70[44],_9B70[54]);
    _9B80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _9B80[i] = -2.16 +0.0025*(i-5);
        }
        //printf("9B %1.2f\n",_9B80[49]);

    _10B40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10B40[i] = -1.92 + 0.0022*(i-5);
            else if (i<90 && i>=45) _10B40[i] = -1.83 + 0.0012*(i-45);
        }
        //printf("10B %1.2f  %1.2f\n",_10B40[44],_10B40[76]);
    _10B50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10B50[i] = -1.95 + 0.0022*(i-5);
            else if (i<90 && i>=45) _10B50[i] = -1.86 + 0.0008*(i-45);
        }
        //printf("10B %1.2f  %1.2f\n",_10B50[44],_10B50[85]);
    _10B60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10B60[i] = -2.01 + 0.0027*(i-5);
            else if (i<90 && i>=45) _10B60[i] = -1.90 + 0.0012*(i-45);
        }
        //printf("10B %1.2f  %1.2f\n",_10B60[44],_10B60[70]);
    _10B70(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _10B70[i] = -2.03 +0.0025*(i-5);
        }
        //printf("10B %1.2f\n",_10B70[49]);
    _10B80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _10B80[i] = -2.08 +0.0032*(i-5);
        }
        //printf("10B %1.2f\n",_10B80[39]);

    _05PB40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _05PB40[i] = -1.87 + 0.0022*(i-5);
            else if (i<90 && i>=45) _05PB40[i] = -1.78 + 0.0015*(i-45);
        }
        //printf("05PB %1.2f  %1.2f\n",_05PB40[44],_05PB40[74]);
    _05PB50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _05PB50[i] = -1.91 + 0.0022*(i-5);
            else if (i<90 && i>=45) _05PB50[i] = -1.82 + 0.001*(i-45);
        }
        //printf("05PB %1.2f  %1.2f\n",_05PB50[44],_05PB50[85]);
    _05PB60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _05PB60[i] = -1.96 + 0.0027*(i-5);
            else if (i<90 && i>=45) _05PB60[i] = -1.85 + 0.0013*(i-45);
        }
        //printf("05PB %1.2f  %1.2f\n",_05PB60[44],_05PB60[70]);
    _05PB70(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _05PB70[i] = -1.99 + 0.0027*(i-5);
            else if (i<90 && i>=45) _05PB70[i] = -1.88 + 0.001*(i-45);
        }
        //printf("05PB %1.2f  %1.2f\n",_05PB70[44],_05PB70[54]);
    _05PB80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _05PB80[i] = -2.03 +0.003*(i-5);
        }
        //printf("05PB %1.2f\n",_05PB80[39]);



    //blue purple correction
    //between 15PB to 4P
    //maximum deviation 75PB

    //15PB
    _15PB10(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _15PB10[i] = -1.66 +0.0035*(i-5);
        }
        //printf("15 %1.2f\n",_15PB10[49]);
    _15PB20(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _15PB20[i] = -1.71 +0.00275*(i-5);
            else if (i<90 && i>=45) _15PB20[i] = -1.60+0.0012*(i-45);
        }
        //printf("15 %1.2f  %1.2f\n",_15PB20[44],_15PB20[89]);

    _15PB30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _15PB30[i] = -1.75 +0.0025*(i-5);
            else if (i<90 && i>=45) _15PB30[i] = -1.65+0.002*(i-45);
        }
        //printf("15 %1.2f  %1.2f\n",_15PB30[44],_15PB30[89]);

    _15PB40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _15PB40[i] = -1.79 +0.002*(i-5);
            else if (i<90 && i>=45) _15PB40[i] = -1.71+0.002*(i-45);
        }
        //printf("15 %1.2f  %1.2f\n",_15PB40[44],_15PB40[89]);

    _15PB50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _15PB50[i] = -1.82 +0.002*(i-5);
            else if (i<90 && i>=45) _15PB50[i] = -1.74+0.0011*(i-45);
        }
        //printf("15 %1.2f  %1.2f\n",_15PB50[44],_15PB50[89]);

    _15PB60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _15PB60[i] = -1.87 +0.0025*(i-5);
            else if (i<90 && i>=45) _15PB60[i] = -1.77+0.001*(i-45);
        }
        //printf("15 %1.2f  %1.2f\n",_15PB60[44],_15PB60[89]);
    _15PB70(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _15PB70[i] = -1.90 +0.0027*(i-5);
        }
        //    printf("15 %1.2f\n",_15PB70[49]);
    _15PB80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _15PB80[i] = -1.93 +0.0027*(i-5);
        }
        //printf("15 %1.2f %1.2f\n",_15PB80[38], _15PB80[49]);

    //3PB
    _3PB10(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _3PB10[i] = -1.56 +0.005*(i-5);
            else if (i<90 && i>=45) _3PB10[i] = -1.36+0.001*(i-45);
        }
        //printf("30 %1.2f  %1.2f\n",_3PB10[44],_3PB10[89]);

    _3PB20(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _3PB20[i] = -1.59 +0.00275*(i-5);
            else if (i<90 && i>=45) _3PB20[i] = -1.48+0.003*(i-45);
        }
        //printf("30 %1.2f  %1.2f\n",_3PB20[44],_3PB20[89]);

    _3PB30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _3PB30[i] = -1.62 +0.00225*(i-5);
            else if (i<90 && i>=45) _3PB30[i] = -1.53+0.0032*(i-45);
        }
        //printf("30 %1.2f  %1.2f\n",_3PB30[44],_3PB30[89]);

    _3PB40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _3PB40[i] = -1.64 +0.0015*(i-5);
            else if (i<90 && i>=45) _3PB40[i] = -1.58+0.0025*(i-45);
        }
        //printf("30 %1.2f  %1.2f\n",_3PB40[44],_3PB40[89]);

    _3PB50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _3PB50[i] = -1.69 +0.00175*(i-5);
            else if (i<90 && i>=45) _3PB50[i] = -1.62+0.002*(i-45);
        }
        //printf("30 %1.2f  %1.2f\n",_3PB50[44],_3PB50[89]);

    _3PB60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _3PB60[i] = -1.73 +0.002*(i-5);
            else if (i<90 && i>=45) _3PB60[i] = -1.65+0.0012*(i-45);
        }
        //printf("30 %1.2f  %1.2f\n",_3PB60[44],_3PB60[89]);
    _3PB70(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _3PB70[i] = -1.76 +0.002*(i-5);
        }
        //printf("30 %1.2f\n",_3PB70[49]);
    _3PB80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _3PB80[i] = -1.78 +0.0025*(i-5);
        }
        //printf("30 %1.2f %1.2f\n",_3PB80[38], _3PB80[49]);

    //45PB
    _45PB10(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _45PB10[i] = -1.46 +0.0045*(i-5);
            else if (i<90 && i>=45) _45PB10[i] = -1.28+0.0025*(i-45);
        }
        //printf("45 %1.2f  %1.2f\n",_45PB10[44],_45PB10[89]);

    _45PB20(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _45PB20[i] = -1.48 +0.00275*(i-5);
            else if (i<90 && i>=45) _45PB20[i] = -1.37+0.0025*(i-45);
        }
        //printf("45 %1.2f  %1.2f\n",_45PB20[44],_45PB20[89]);

    _45PB30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _45PB30[i] = -1.51 +0.00175*(i-5);
            else if (i<90 && i>=45) _45PB30[i] = -1.44+0.0035*(i-45);
        }
        //printf("45 %1.2f  %1.2f\n",_45PB30[44],_45PB30[89]);

    _45PB40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _45PB40[i] = -1.52 +0.001*(i-5);
            else if (i<90 && i>=45) _45PB40[i] = -1.48+0.003*(i-45);
        }
        //printf("45 %1.2f  %1.2f\n",_45PB40[44],_45PB40[89]);

    _45PB50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _45PB50[i] = -1.55 +0.001*(i-5);
            else if (i<90 && i>=45) _45PB50[i] = -1.51+0.0022*(i-45);
        }
        //printf("45 %1.2f  %1.2f\n",_45PB50[44],_45PB50[89]);

    _45PB60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _45PB60[i] = -1.6 +0.0015*(i-5);
            else if (i<90 && i>=45) _45PB60[i] = -1.54+0.001*(i-45);
        }
        //printf("45 %1.2f  %1.2f\n",_45PB60[44],_45PB60[89]);
    _45PB70(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _45PB70[i] = -1.63 +0.0017*(i-5);
        }
        //printf("45 %1.2f\n",_45PB70[49]);
    _45PB80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _45PB80[i] = -1.67 +0.0025*(i-5);
        }
        //printf("45 %1.2f %1.2f\n",_45PB80[38], _45PB80[49]);

    //_6PB
    _6PB10(maxInd);
        for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
            if (i<45 && i>5) _6PB10[i] = -1.33 +0.005*(i-5);
            else if (i<85 && i>=45) _6PB10[i] = -1.13+0.0045*(i-45);
            else if (i<140 && i >=85) _6PB10[i] = -0.95+0.0015*(i-85);
        }
        //printf("60 %1.2f  %1.2f %1.2f\n",_6PB10[44],_6PB10[84],_6PB10[139]);

    _6PB20(maxInd);
        for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
            if (i<45 && i>5) _6PB20[i] = -1.36 +0.004*(i-5);
            else if (i<85 && i>=45) _6PB20[i] = -1.20+0.00375*(i-45);
            else if (i<140 && i >=85) _6PB20[i] = -1.05+0.0017*(i-85);
        }
        //printf("60 %1.2f  %1.2f %1.2f\n",_6PB20[44],_6PB20[84],_6PB20[139]);

    _6PB30(maxInd);
        for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
            if (i<45 && i>5) _6PB30[i] = -1.38 +0.00225*(i-5);
            else if (i<85 && i>=45) _6PB30[i] = -1.29+0.00375*(i-45);
            else if (i<140 && i >=85) _6PB30[i] = -1.14+0.002*(i-85);
        }
        //printf("60 %1.2f  %1.2f %1.2f\n",_6PB30[44],_6PB30[84],_6PB30[139]);

    _6PB40(maxInd);
        for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
            if (i<45 && i>5) _6PB40[i] = -1.39 +0.00125*(i-5);
            else if (i<85 && i>=45) _6PB40[i] = -1.34+0.00275*(i-45);
            else if (i<140 && i >=85) _6PB40[i] = -1.23+0.002*(i-85);
        }
        //printf("60 %1.2f  %1.2f %1.2f\n",_6PB40[44],_6PB40[84],_6PB40[139]);

    _6PB50(maxInd2);//limits  -1.3   -1.11
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _6PB50[i] = -1.43 +0.00125*(i-5);
            else if (i<90 && i>=45) _6PB50[i] = -1.38+0.00225*(i-45);
        }
        //printf("60 %1.2f  %1.2f \n",_6PB50[44],_6PB50[89]);

    _6PB60(maxInd2);//limits  -1.3   -1.11
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _6PB60[i] = -1.46 +0.0012*(i-5);
            else if (i<90 && i>=45) _6PB60[i] = -1.40+0.000875*(i-45);
        }
        //printf("60 %1.2f  %1.2f\n",_6PB60[44],_6PB60[89]);
    _6PB70(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _6PB70[i] = -1.49 +0.0018*(i-5);
        }
        //printf("6 %1.2f\n",_6PB70[49]);
    _6PB80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _6PB80[i] = -1.52 +0.0022*(i-5);
        }
        //printf("6 %1.2f %1.2f\n",_6PB80[38], _6PB80[49]);


    //_75PB : notation Munsell for maximum deviation blue purple
    _75PB10(maxInd);//limits hue -1.23  -0.71  _75PBx   x=Luminance  eg_75PB10 for L >5 and L<=15
        for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
            if (i<45 && i>5) _75PB10[i] = -1.23 +0.0065*(i-5);
            else if (i<85 && i>=45) _75PB10[i] = -0.97+0.00375*(i-45);
            else if (i<140 && i >=85) _75PB10[i] = -0.82+0.0015*(i-85);
        }
        //printf("75 %1.2f  %1.2f %1.2f\n",_75PB10[44],_75PB10[84],_75PB10[139]);

    _75PB20(maxInd);//limits -1.24  -0.79  for L>15 <=25
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _75PB20[i] = -1.24 +0.004*(i-5);
            else if (i<85 && i>=45) _75PB20[i] = -1.08+0.00425*(i-45);
            else if (i<140 && i >=85) _75PB20[i] = -0.91+0.0017*(i-85);
        }
        //printf("75 %1.2f  %1.2f %1.2f\n",_75PB20[44],_75PB20[84],_75PB20[139]);

    _75PB30(maxInd);//limits -1.25  -0.85
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _75PB30[i] = -1.25 +0.00275*(i-5);
            else if (i<85 && i>=45) _75PB30[i] = -1.14+0.004*(i-45);
            else if (i<140 && i >=85) _75PB30[i] = -0.98+0.0015*(i-85);
        }
        //printf("75 %1.2f  %1.2f %1.2f\n",_75PB30[44],_75PB30[84],_75PB30[139]);

    _75PB40(maxInd);//limits  -1.27  -0.92
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _75PB40[i] = -1.27 +0.002*(i-5);
            else if (i<85 && i>=45) _75PB40[i] = -1.19+0.003*(i-45);
            else if (i<140 && i >=85) _75PB40[i] = -1.07+0.0022*(i-85);
        }
        //printf("75 %1.2f  %1.2f %1.2f\n",_75PB40[44],_75PB40[84],_75PB40[139]);

    _75PB50(maxInd2);//limits  -1.3   -1.11
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _75PB50[i] = -1.3 +0.00175*(i-5);
            else if (i<90 && i>=45) _75PB50[i] = -1.23+0.0025*(i-45);
        }
        //printf("75 %1.2f  %1.2f\n",_75PB50[44],_75PB50[89]);

    _75PB60(maxInd2);
        for (int i=0; i<maxInd2; i++) {//limits -1.32  -1.17
            if (i<45 && i>5) _75PB60[i] = -1.32 +0.0015*(i-5);
            else if (i<90 && i>=45) _75PB60[i] = -1.26+0.002*(i-45);
        }
        //printf("75 %1.2f  %1.2f \n",_75PB60[44],_75PB60[89]);

    _75PB70(maxInd3);
        for (int i=0; i<maxInd3; i++) {//limits  -1.34  -1.27
            if (i<50 && i>5) _75PB70[i] = -1.34 +0.002*(i-5);
        }
    _75PB80(maxInd3);
        for (int i=0; i<maxInd3; i++) {//limits -1.35  -1.29
            if (i<50 && i>5) _75PB80[i] = -1.35 +0.00125*(i-5);
        }


    _9PB10(maxInd);
        for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
            if (i<45 && i>5) _9PB10[i] = -1.09 +0.00475*(i-5);
            else if (i<85 && i>=45) _9PB10[i] = -0.9+0.003*(i-45);
            else if (i<140 && i >=85) _9PB10[i] = -0.78+0.0013*(i-85);
        }
        //printf("90 %1.2f  %1.2f %1.2f\n",_9PB10[44],_9PB10[84],_9PB10[139]);

    _9PB20(maxInd);
        for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
            if (i<45 && i>5) _9PB20[i] = -1.12 +0.0035*(i-5);
            else if (i<85 && i>=45) _9PB20[i] = -0.98+0.00325*(i-45);
            else if (i<140 && i >=85) _9PB20[i] = -0.85+0.0015*(i-85);
        }
        //printf("90 %1.2f  %1.2f %1.2f\n",_9PB20[44],_9PB20[84],_9PB20[139]);

    _9PB30(maxInd);
        for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
            if (i<45 && i>5) _9PB30[i] = -1.14 +0.0028*(i-5);
            else if (i<85 && i>=45) _9PB30[i] = -1.03+0.003*(i-45);
            else if (i<140 && i >=85) _9PB30[i] = -0.91+0.0017*(i-85);
        }
        //printf("90 %1.2f  %1.2f %1.2f\n",_9PB30[44],_9PB30[84],_9PB30[139]);

    _9PB40(maxInd);
        for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
            if (i<45 && i>5) _9PB40[i] = -1.16 +0.002*(i-5);
            else if (i<85 && i>=45) _9PB40[i] = -1.08+0.00275*(i-45);
            else if (i<140 && i >=85) _9PB40[i] = -0.97+0.0016*(i-85);
        }
        //printf("90 %1.2f  %1.2f %1.2f\n",_9PB40[44],_9PB40[84],_9PB40[139]);

    _9PB50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _9PB50[i] = -1.19 +0.00175*(i-5);
            else if (i<90 && i>=45) _9PB50[i] = -1.12+0.00225*(i-45);
        }
        //printf("90 %1.2f  %1.2f \n",_9PB50[44],_9PB50[84]);

    _9PB60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _9PB60[i] = -1.21 +0.0015*(i-5);
            else if (i<90 && i>=45) _9PB60[i] = -1.15+0.002*(i-45);
        }
        //printf("90 %1.2f  %1.2f \n",_9PB60[44],_9PB60[89]);
    _9PB70(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _9PB70[i] = -1.23 +0.0018*(i-5);
        }
            //printf("9 %1.2f\n",_9PB70[49]);
    _9PB80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _9PB80[i] = -1.24 +0.002*(i-5);
        }
            //printf("9 %1.2f %1.2f\n",_9PB80[38], _9PB80[49]);


    //10PB
    _10PB10(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10PB10[i] = -1.02 +0.00425*(i-5);
            else if (i<85 && i>=45) _10PB10[i] = -0.85+0.0025*(i-45);
            else if (i<140 && i >=85) _10PB10[i] = -0.75+0.0012*(i-85);
        }
        //printf("10 %1.2f  %1.2f %1.2f\n",_10PB10[44],_10PB10[84],_10PB10[139]);

    _10PB20(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10PB20[i] = -1.05 +0.00325*(i-5);
            else if (i<85 && i>=45) _10PB20[i] = -0.92+0.00275*(i-45);
            else if (i<140 && i >=85) _10PB20[i] = -0.81+0.0014*(i-85);
        }
        //printf("10 %1.2f  %1.2f %1.2f\n",_10PB20[44],_10PB20[84],_10PB20[139]);

    _10PB30(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10PB30[i] = -1.07 +0.00275*(i-5);
            else if (i<85 && i>=45) _10PB30[i] = -0.96+0.0025*(i-45);
            else if (i<140 && i >=85) _10PB30[i] = -0.86+0.0015*(i-85);
        }
        //printf("10 %1.2f  %1.2f %1.2f\n",_10PB30[44],_10PB30[84],_10PB30[139]);

    _10PB40(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10PB40[i] = -1.09 +0.002*(i-5);
            else if (i<85 && i>=45) _10PB40[i] = -1.01+0.00225*(i-45);
            else if (i<140 && i >=85) _10PB40[i] = -0.92+0.0016*(i-85);
        }
        //printf("10 %1.2f  %1.2f %1.2f\n",_10PB40[44],_10PB40[84],_10PB40[139]);

    _10PB50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10PB50[i] = -1.12 +0.00175*(i-5);
            else if (i<90 && i>=45) _10PB50[i] = -1.05+0.00225*(i-45);
        }
        //printf("10 %1.2f  %1.2f\n",_10PB50[44],_10PB50[84]);

    _10PB60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10PB60[i] = -1.14 +0.0015*(i-5);
            else if (i<90 && i>=45) _10PB60[i] = -1.08+0.00225*(i-45);
        }
        //printf("10 %1.2f  %1.2f\n",_10PB60[44],_10PB60[89]);


    //1P
    _1P10(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1P10[i] = -0.96 +0.00375*(i-5);
            else if (i<85 && i>=45) _1P10[i] = -0.81+0.00225*(i-45);
            else if (i<140 && i >=85) _1P10[i] = -0.72+0.001*(i-85);
        }
        //printf("1P %1.2f  %1.2f %1.2f\n",_1P10[44],_1P10[84],_1P10[139]);

    _1P20(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1P20[i] = -1.0 +0.00325*(i-5);
            else if (i<85 && i>=45) _1P20[i] = -0.87+0.0025*(i-45);
            else if (i<140 && i >=85) _1P20[i] = -0.77+0.0012*(i-85);
        }
        //printf("1P %1.2f  %1.2f %1.2f\n",_1P20[44],_1P20[84],_1P20[139]);

    _1P30(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1P30[i] = -1.02 +0.00275*(i-5);
            else if (i<85 && i>=45) _1P30[i] = -0.91+0.00225*(i-45);
            else if (i<140 && i >=85) _1P30[i] = -0.82+0.0011*(i-85);
        }
        //printf("1P %1.2f  %1.2f %1.2f\n",_1P30[44],_1P30[84],_1P30[139]);

    _1P40(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1P40[i] = -1.04 +0.00225*(i-5);
            else if (i<85 && i>=45) _1P40[i] = -0.95+0.00225*(i-45);
            else if (i<140 && i >=85) _1P40[i] = -0.86+0.0015*(i-85);
        }
        //printf("1P %1.2f  %1.2f %1.2f\n",_1P40[44],_1P40[84],_1P40[139]);

    _1P50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _1P50[i] = -1.06 +0.002*(i-5);
            else if (i<90 && i>=45) _1P50[i] = -0.98+0.00175*(i-45);
        }
        //printf("1P %1.2f  %1.2f \n",_1P50[44],_1P50[89]);

    _1P60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _1P60[i] = -1.07 +0.0015*(i-5);
            else if (i<90 && i>=45) _1P60[i] = -1.01+0.00175*(i-45);
        }
        //printf("1P %1.2f  %1.2f \n",_1P60[44],_1P60[84],_1P60[139]);

            //4P
    _4P10(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _4P10[i] = -0.78 +0.002*(i-5);
            else if (i<85 && i>=45) _4P10[i] = -0.7+0.00125*(i-45);
            else if (i<140 && i >=85) _4P10[i] = -0.65+0.001*(i-85);
        }
        //printf("4P %1.2f  %1.2f %1.2f\n",_4P10[44],_4P10[84],_4P10[139]);

    _4P20(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _4P20[i] = -0.84 +0.0025*(i-5);
            else if (i<85 && i>=45) _4P20[i] = -0.74+0.00175*(i-45);
            else if (i<140 && i >=85) _4P20[i] = -0.67+0.00085*(i-85);
        }
        //printf("4P %1.2f  %1.2f %1.2f\n",_4P20[44],_4P20[84],_4P20[139]);

    _4P30(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _4P30[i] = -0.85 +0.00225*(i-5);
            else if (i<85 && i>=45) _4P30[i] = -0.76+0.00125*(i-45);
            else if (i<140 && i >=85) _4P30[i] = -0.71+0.001*(i-85);
        }
        //printf("4P %1.2f  %1.2f %1.2f\n",_4P30[44],_4P30[84],_4P30[139]);

    _4P40(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _4P40[i] = -0.87 +0.00175*(i-5);
            else if (i<85 && i>=45) _4P40[i] = -0.8+0.00175*(i-45);
            else if (i<140 && i >=85) _4P40[i] = -0.73+0.00075*(i-85);
        }
        //printf("4P %1.2f  %1.2f %1.2f\n",_4P40[44],_4P40[84],_4P40[139]);

    _4P50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _4P50[i] = -0.88 +0.0015*(i-5);
            else if (i<90 && i>=45) _4P50[i] = -0.82+0.0015*(i-45);
        }
        //printf("4P %1.2f  %1.2f \n",_4P50[44],_4P50[89]);

    _4P60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _4P60[i] = -0.89 +0.00125*(i-5);
            else if (i<90 && i>=45) _4P60[i] = -0.84+0.00125*(i-45);
        }
        //printf("4P %1.2f  %1.2f\n",_4P60[44],_4P60[89]);


    //red yellow correction
    _10YR20(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10YR20[i] = 1.22 +0.002*(i-5);
            else if (i<90 && i>=45) _10YR20[i] = 1.30+0.006*(i-45);
        }
        //printf("10YR  %1.2f  %1.2f\n",_10YR20[44],_10YR20[56]);
    _10YR30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10YR30[i] = 1.27 +0.00175*(i-5);
            else if (i<90 && i>=45) _10YR30[i] = 1.34+0.0017*(i-45);
        }
        //printf("10YR  %1.2f  %1.2f\n",_10YR30[44],_10YR30[75]);
    _10YR40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10YR40[i] = 1.32 +0.00025*(i-5);
            else if (i<90 && i>=45) _10YR40[i] = 1.33+0.0015*(i-45);
        }
        //printf("10YR  %1.2f  %1.2f\n",_10YR40[44],_10YR40[85]);
    _10YR50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10YR50[i] = 1.35 +0.000*(i-5);
            else if (i<90 && i>=45) _10YR50[i] = 1.35+0.0012*(i-45);
        }
        //printf("10YR  %1.2f  %1.2f\n",_10YR50[44],_10YR50[80]);
    _10YR60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10YR60[i] = 1.38 - 0.00025*(i-5);
            else if (i<85 && i>=45) _10YR60[i] = 1.37+0.0005*(i-45);
            else if (i<140 && i >=85) _10YR60[i] = 1.39+0.0013*(i-85);
        }
        //printf("10YR  %1.2f  %1.2f %1.2f\n",_10YR60[44],_10YR60[85],_10YR60[139] );
    _10YR70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10YR70[i] = 1.41 - 0.0005*(i-5);
            else if (i<85 && i>=45) _10YR70[i] = 1.39+0.000*(i-45);
            else if (i<140 && i >=85) _10YR70[i] = 1.39+0.0013*(i-85);
        }
        //printf("10YR  %1.2f  %1.2f %1.2f\n",_10YR70[44],_10YR70[85],_10YR70[139] );
    _10YR80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10YR80[i] = 1.45 - 0.00125*(i-5);
            else if (i<85 && i>=45) _10YR80[i] = 1.40+0.000*(i-45);
            else if (i<140 && i >=85) _10YR80[i] = 1.40+0.00072*(i-85);//1.436
        }
        //printf("10YR  %1.2f  %1.2f %1.2f\n",_10YR80[44],_10YR80[84],_10YR80[139] );
    _10YR90(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10YR90[i] = 1.48 -0.001*(i-5);
            else if (i<90 && i>=45) _10YR90[i] = 1.44-0.0009*(i-45);
        }
        //printf("10YR  %1.2f  %1.2f\n",_10YR90[45],_10YR90[80]);
    _85YR20(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _85YR20[i] = 1.12 +0.004*(i-5);
        }

        //printf("85YR  %1.2f \n",_85YR20[44]);
    _85YR30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _85YR30[i] = 1.16 + 0.0025*(i-5);
            else if (i<90 && i>=45) _85YR30[i] = 1.26+0.0028*(i-45);
        }
        //printf("85YR  %1.2f  %1.2f\n",_85YR30[44],_85YR30[75]);
    _85YR40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _85YR40[i] = 1.20 + 0.0015*(i-5);
            else if (i<90 && i>=45) _85YR40[i] = 1.26+0.0024*(i-45);
        }
        //printf("85YR  %1.2f  %1.2f\n",_85YR40[44],_85YR40[75]);
    _85YR50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _85YR50[i] = 1.24 + 0.0005*(i-5);
            else if (i<85 && i>=45) _85YR50[i] = 1.26+0.002*(i-45);
            else if (i<140 && i >=85) _85YR50[i] = 1.34+0.0015*(i-85);
        }
        //printf("85YR  %1.2f  %1.2f %1.2f\n",_85YR50[44],_85YR50[85],_85YR50[110] );
    _85YR60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _85YR60[i] = 1.27 + 0.00025*(i-5);
            else if (i<85 && i>=45) _85YR60[i] = 1.28+0.0015*(i-45);
            else if (i<140 && i >=85) _85YR60[i] = 1.34+0.0012*(i-85);
        }
        //printf("85YR  %1.2f  %1.2f %1.2f\n",_85YR60[44],_85YR60[85],_85YR60[139] );

    _85YR70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _85YR70[i] = 1.31 - 0.00025*(i-5);
            else if (i<85 && i>=45) _85YR70[i] = 1.30+0.0005*(i-45);
            else if (i<140 && i >=85) _85YR70[i] = 1.32+0.0012*(i-85);
        }
        //printf("85YR  %1.2f  %1.2f %1.2f\n",_85YR70[44],_85YR70[85],_85YR70[139] );
    _85YR80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _85YR80[i] = 1.35 - 0.00075*(i-5);
            else if (i<85 && i>=45) _85YR80[i] = 1.32+0.00025*(i-45);
            else if (i<140 && i >=85) _85YR80[i] = 1.33+0.00125*(i-85);
        }
        //printf("85YR  %1.2f  %1.2f %1.2f\n",_85YR80[44],_85YR80[85],_85YR80[139] );
    _85YR90(maxInd2);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _85YR90[i] = 1.39 - 0.00125*(i-5);
            else if (i<90 && i>=45) _85YR90[i] = 1.34+0.00*(i-45);
        }
        //printf("85YR  %1.2f  %1.2f\n",_85YR90[44],_85YR90[85]);

    //7YR
    _7YR30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7YR30[i] = 1.06 + 0.0028*(i-5);
            else if (i<90 && i>=45) _7YR30[i] = 1.17+0.0045*(i-45);
        }
        //printf("7YR  %1.2f  %1.2f\n",_7YR30[44],_7YR30[66]);
    _7YR40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7YR40[i] = 1.10 + 0.0018*(i-5);
            else if (i<90 && i>=45) _7YR40[i] = 1.17+0.0035*(i-45);
        }
        //printf("7YR  %1.2f  %1.2f\n",_7YR40[44],_7YR40[89]);
    _7YR50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7YR50[i] = 1.14 + 0.00125*(i-5);
            else if (i<90 && i>=45) _7YR50[i] = 1.19+0.002*(i-45);
        }
        //printf("7YR  %1.2f  %1.2f\n",_7YR50[44],_7YR50[89] );
    _7YR60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7YR60[i] = 1.17 + 0.00075*(i-5);
            else if (i<85 && i>=45) _7YR60[i] = 1.20+0.00175*(i-45);
            else if (i<140 && i >=85) _7YR60[i] = 1.27+0.002*(i-85);
        }
        //printf("7YR  %1.2f  %1.2f %1.2f\n",_7YR60[44],_7YR60[84],_7YR60[125] );

    _7YR70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7YR70[i] = 1.20 + 0.0005*(i-5);
            else if (i<85 && i>=45) _7YR70[i] = 1.22+0.00125*(i-45);
            else if (i<140 && i >=85) _7YR70[i] = 1.27+0.0015*(i-85);
        }
        //printf("7YR  %1.2f  %1.2f %1.2f\n",_7YR70[44],_7YR70[84],_7YR70[125] );
    _7YR80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _7YR80[i] = 1.29 - 0.0008*(i-5);
        }
        //printf("7YR  %1.2f \n",_7YR80[44] );
    _55YR30(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _55YR30[i] = 0.96 + 0.0038*(i-5);
        }
        //printf("55YR  %1.2f \n",_55YR30[44] );
    _55YR40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _55YR40[i] = 1.01 + 0.0022*(i-5);
            else if (i<90 && i>=45) _55YR40[i] = 1.10+0.0037*(i-45);
        }
        //printf("55YR  %1.2f  %1.2f\n",_55YR40[44],_55YR40[89] );
    _55YR50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _55YR50[i] = 1.06 + 0.0015*(i-5);
            else if (i<85 && i>=45) _55YR50[i] = 1.12+0.00225*(i-45);
            else if (i<140 && i >=85) _55YR50[i] = 1.21+0.0015*(i-85);
        }
        //printf("55YR  %1.2f  %1.2f %1.2f\n",_55YR50[44],_55YR50[84],_55YR50[125] );
    _55YR60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _55YR60[i] = 1.08 + 0.0012*(i-5);
            else if (i<85 && i>=45) _55YR60[i] = 1.13+0.0018*(i-45);
            else if (i<140 && i >=85) _55YR60[i] = 1.20+0.0025*(i-85);
        }
        //printf("55YR  %1.2f  %1.2f %1.2f\n",_55YR60[44],_55YR60[84],_55YR60[125] );
    _55YR70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _55YR70[i] = 1.11 + 0.00075*(i-5);
            else if (i<85 && i>=45) _55YR70[i] = 1.14+0.0012*(i-45);
            else if (i<140 && i >=85) _55YR70[i] = 1.19+0.00225*(i-85);
        }
        //printf("55YR  %1.2f  %1.2f %1.2f\n",_55YR70[44],_55YR70[84],_55YR70[125] );
    _55YR80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _55YR80[i] = 1.16 + 0.00*(i-5);
            else if (i<85 && i>=45) _55YR80[i] = 1.16+0.00075*(i-45);
            else if (i<140 && i >=85) _55YR80[i] = 1.19+0.00175*(i-85);
        }
        //printf("55YR  %1.2f  %1.2f %1.2f\n",_55YR80[44],_55YR80[84],_55YR80[125] );
    _55YR90(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _55YR90[i] = 1.19 - 0.0005*(i-5);
        }
        //printf("55YR  %1.2f \n",_55YR90[44] );

    _4YR30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _4YR30[i] = 0.87 + 0.0035*(i-5);
            else if (i<90 && i>=45) _4YR30[i] = 1.01+0.0043*(i-45);
        }
        //printf("4YR  %1.2f  %1.2f\n",_4YR30[44],_4YR30[78] );
    _4YR40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _4YR40[i] = 0.92 + 0.0025*(i-5);
            else if (i<90 && i>=45) _4YR40[i] = 1.02+0.0033*(i-45);
        }
        //printf("4YR  %1.2f  %1.2f\n",_4YR40[44],_4YR40[74] );
    _4YR50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _4YR50[i] = 0.97 + 0.0015*(i-5);
            else if (i<90 && i>=45) _4YR50[i] = 1.03+0.0025*(i-45);
        }
        //printf("4YR  %1.2f  %1.2f\n",_4YR50[44],_4YR50[85] );
    _4YR60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _4YR60[i] = 0.99 + 0.00125*(i-5);
            else if (i<85 && i>=45) _4YR60[i] = 1.04+0.002*(i-45);
            else if (i<140 && i >=85) _4YR60[i] = 1.12+0.003*(i-85);
        }
        //printf("4YR  %1.2f  %1.2f %1.2f\n",_4YR60[44],_4YR60[84],_4YR60[125] );
    _4YR70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _4YR70[i] = 1.02 + 0.00075*(i-5);
            else if (i<85 && i>=45) _4YR70[i] = 1.05+0.00175*(i-45);
            else if (i<140 && i >=85) _4YR70[i] = 1.12+0.002*(i-85);
        }
        //printf("4YR  %1.2f  %1.2f %1.2f\n",_4YR70[44],_4YR70[84],_4YR70[125] );
    _4YR80(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<50 && i>5) _4YR80[i] = 1.09 - 0.0002*(i-5);
        }
        //printf("4YR  %1.2f \n",_4YR80[41] );

    _25YR30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _25YR30[i] = 0.77 + 0.004*(i-5);
            else if (i<90 && i>=45) _25YR30[i] = 0.94+0.004*(i-45);
        }
        //printf("25YR  %1.2f  %1.2f\n",_25YR30[44],_25YR30[74] );
    _25YR40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _25YR40[i] = 0.82 + 0.003*(i-5);
            else if (i<90 && i>=45) _25YR40[i] = 0.94+0.002*(i-45);
        }
        //printf("25YR  %1.2f  %1.2f\n",_25YR40[44],_25YR40[84] );
    _25YR50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _25YR50[i] = 0.87+ 0.002*(i-5);
            else if (i<90 && i>=45) _25YR50[i] = 0.95+0.003*(i-45);
        }
        //printf("25YR  %1.2f  %1.2f\n",_25YR50[44],_25YR50[84] );
    _25YR60(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _25YR60[i] = 0.89+ 0.0015*(i-5);
            else if (i<90 && i>=45) _25YR60[i] = 0.95+0.004*(i-45);
        }
        //printf("25YR  %1.2f  %1.2f\n",_25YR60[44],_25YR60[84] );
    _25YR70(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _25YR70[i] = 0.92+ 0.001*(i-5);
            else if (i<90 && i>=45) _25YR70[i] = 0.96+0.003*(i-45);
        }
        //printf("25YR  %1.2f  %1.2f\n",_25YR70[44],_25YR70[84] );

    _10R30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10R30[i] = 0.62 + 0.00225*(i-5);
            else if (i<90 && i>=45) _10R30[i] = 0.71+0.003*(i-45);
        }
        //printf("10R  %1.2f  %1.2f\n",_10R30[44],_10R30[84] );
    _10R40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10R40[i] = 0.66 + 0.0025*(i-5);
            else if (i<90 && i>=45) _10R40[i] = 0.76+0.0035*(i-45);
        }
        //printf("10R  %1.2f  %1.2f\n",_10R40[44],_10R40[84] );
    _10R50(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10R50[i] = 0.71 + 0.002*(i-5);
            else if (i<90 && i>=45) _10R50[i] = 0.79+0.0043*(i-45);
        }
        //printf("10R  %1.2f  %1.2f\n",_10R50[44],_10R50[84] );
    _10R60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10R60[i] = 0.73 + 0.00175*(i-5);
            else if (i<85 && i>=45) _10R60[i] = 0.80 +0.0033*(i-45);
            else if (i<140 && i >=85) _10R60[i] = 0.93+0.0018*(i-85);
        }
        //printf("10R  %1.2f  %1.2f %1.2f\n",_10R60[44],_10R60[84],_10R60[125] );
    _10R70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10R70[i] = 0.75 + 0.0015*(i-5);
            else if (i<85 && i>=45) _10R70[i] = 0.81 +0.0017*(i-45);
            else if (i<140 && i >=85) _10R70[i] = 0.88+0.0025*(i-85);
        }
        //printf("10R  %1.2f  %1.2f %1.2f\n",_10R70[44],_10R70[84],_10R70[125] );

    _9R30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _9R30[i] = 0.57 + 0.002*(i-5);
            else if (i<90 && i>=45) _9R30[i] = 0.65+0.0018*(i-45);
        }
        //printf("9R  %1.2f  %1.2f\n",_9R30[44],_9R30[84] );
    _9R40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _9R40[i] = 0.61 + 0.002*(i-5);
            else if (i<90 && i>=45) _9R40[i] = 0.69+0.0025*(i-45);
        }
        //printf("9R  %1.2f  %1.2f\n",_9R40[44],_9R40[84] );
    _9R50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _9R50[i] = 0.66 + 0.00175*(i-5);
            else if (i<85 && i>=45) _9R50[i] = 0.73 +0.0025*(i-45);
            else if (i<140 && i >=85) _9R50[i] = 0.83+0.0035*(i-85);
        }
        //printf("9R  %1.2f  %1.2f %1.2f\n",_9R50[44],_9R50[84],_9R50[125] );
    _9R60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _9R60[i] = 0.68 + 0.0015*(i-5);
            else if (i<85 && i>=45) _9R60[i] = 0.74 +0.0022*(i-45);
            else if (i<140 && i >=85) _9R60[i] = 0.93+0.0022*(i-85);
        }
        //printf("9R  %1.2f  %1.2f %1.2f\n",_9R60[44],_9R60[84],_9R60[125] );
    _9R70(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _9R70[i] = 0.70 + 0.0012*(i-5);
            else if (i<90 && i>=45) _9R70[i] = 0.75+0.0013*(i-45);
        }
        //printf("9R  %1.2f  %1.2f\n",_9R70[44],_9R70[84] );

    _7R30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7R30[i] = 0.48 + 0.0015*(i-5);
            else if (i<90 && i>=45) _7R30[i] = 0.54-0.0005*(i-45);
        }
        //printf("7R  %1.2f  %1.2f\n",_7R30[44],_7R30[84] );
    _7R40(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7R40[i] = 0.51 + 0.0015*(i-5);
            else if (i<90 && i>=45) _7R40[i] = 0.57+0.0005*(i-45);
        }
        //printf("7R  %1.2f  %1.2f\n",_7R40[44],_7R40[84] );
    _7R50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7R50[i] = 0.54 + 0.0015*(i-5);
            else if (i<85 && i>=45) _7R50[i] = 0.60 +0.0005*(i-45);
            else if (i<140 && i >=85) _7R50[i] = 0.62+0.0025*(i-85);
        }
        //printf("7R  %1.2f  %1.2f %1.2f\n",_7R50[44],_7R50[84],_7R50[125] );
    _7R60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7R60[i] = 0.58 + 0.00075*(i-5);
            else if (i<85 && i>=45) _7R60[i] = 0.61 +0.00075*(i-45);
            else if (i<140 && i >=85) _7R60[i] = 0.64+0.001*(i-85);
        }
        //printf("7R  %1.2f  %1.2f %1.2f\n",_7R60[44],_7R60[84],_7R60[107] );
    _7R70(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _7R70[i] = 0.59 + 0.00075*(i-5);
            else if (i<90 && i>=45) _7R70[i] = 0.62+0.00075*(i-45);
        }
        //printf("7R  %1.2f  %1.2f\n",_7R70[44],_7R70[84] );

    //5R 1 2 3

    //5R
    _5R10(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _5R10[i] = 0.10 - 0.0018*(i-5);
            else if (i<90 && i>=45) _5R10[i] = 0.035-0.003*(i-45);
        }
        //printf("5R  %1.2f  %1.2f\n",_5R10[44],_5R10[51] );
    _5R20(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _5R20[i] = 0.26 - 0.00075*(i-5);
            else if (i<90 && i>=45) _5R20[i] = 0.023-0.0002*(i-45);
        }
        //printf("5R  %1.2f  %1.2f\n",_5R20[44],_5R20[70] );
    _5R30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _5R30[i] = 0.39 + 0.00075*(i-5);
            else if (i<90 && i>=45) _5R30[i] = 0.42-0.0007*(i-45);
        }
        //printf("5R  %1.2f  %1.2f\n",_5R30[44],_5R30[85] );

    //25R
    _25R10(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<45 && i>5) _25R10[i] = -0.03 - 0.002*(i-5);
        }
        //printf("25R  %1.2f \n",_25R10[44]);
    _25R20(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _25R20[i] = 0.13 - 0.0012*(i-5);
            else if (i<90 && i>=45) _25R20[i] = 0.08-0.002*(i-45);
        }
        //printf("25R  %1.2f  %1.2f\n",_25R20[44],_25R20[69] );
        //25R30: 0.28, 0.26, 0.22
    _25R30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _25R30[i] = 0.28 - 0.0005*(i-5);
            else if (i<90 && i>=45) _25R30[i] = 0.26-0.0009*(i-45);
        }
        //printf("25R  %1.2f  %1.2f\n",_25R30[44],_25R30[85] );


    _10RP10(maxInd3);
        for (int i=0; i<maxInd3; i++) {
            if (i<45 && i>5) _10RP10[i] = -0.16 - 0.0017*(i-5);
        }
        //printf("10RP  %1.2f \n",_10RP10[44]);
    _10RP20(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10RP20[i] = 0.0 - 0.0018*(i-5);
            else if (i<90 && i>=45) _10RP20[i] = -0.07-0.0012*(i-45);
        }
        //printf("10RP  %1.2f  %1.2f\n",_10RP20[44],_10RP20[69] );
    _10RP30(maxInd2);
        for (int i=0; i<maxInd2; i++) {
            if (i<45 && i>5) _10RP30[i] = 0.15 - 0.001*(i-5);
            else if (i<90 && i>=45) _10RP30[i] = 0.11-0.0012*(i-45);
        }
        //printf("10RP  %1.2f  %1.2f\n",_10RP30[44],_10RP30[85] );

    //7G
    _7G30(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7G30[i] = 2.90 + 0.0027*(i-5);
            else if (i<85 && i>=45) _7G30[i] = 3.01+0.0005*(i-45);
            else if (i<140 && i >=85) _7G30[i] = 3.03+0.00075*(i-85);
        }
        //printf("7G  %1.2f  %1.2f %1.2f\n",_7G30[44],_7G30[84],_7G30[125] );
    _7G40(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7G40[i] = 2.89 + 0.00125*(i-5);
            else if (i<85 && i>=45) _7G40[i] = 2.94+0.0015*(i-45);
            else if (i<140 && i >=85) _7G40[i] = 3.0+0.001*(i-85);
        }
        //printf("7G  %1.2f  %1.2f %1.2f\n",_7G40[44],_7G40[84],_7G40[125] );
    _7G50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7G50[i] = 2.87 + 0.0015*(i-5);
            else if (i<85 && i>=45) _7G50[i] = 2.93+0.00125*(i-45);
            else if (i<140 && i >=85) _7G50[i] = 2.98+0.001*(i-85);
        }
        //printf("7G  %1.2f  %1.2f %1.2f\n",_7G50[44],_7G50[84],_7G50[125] );
    _7G60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7G60[i] = 2.86 + 0.00125*(i-5);
            else if (i<85 && i>=45) _7G60[i] = 2.91+0.00125*(i-45);
            else if (i<140 && i >=85) _7G60[i] = 2.96+0.00075*(i-85);
        }
        //printf("7G  %1.2f  %1.2f %1.2f\n",_7G60[44],_7G60[84],_7G60[125] );
    _7G70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7G70[i] = 2.85 + 0.001*(i-5);
            else if (i<85 && i>=45) _7G70[i] = 2.89+0.00125*(i-45);
            else if (i<140 && i >=85) _7G70[i] = 2.94+0.00075*(i-85);
        }
        //printf("7G  %1.2f  %1.2f %1.2f\n",_7G70[44],_7G70[84],_7G70[125] );
    _7G80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _7G80[i] = 2.84 + 0.001*(i-5);
            else if (i<85 && i>=45) _7G80[i] = 2.88+0.001*(i-45);
            else if (i<140 && i >=85) _7G80[i] = 2.92+0.001*(i-85);
        }
        //printf("7G  %1.2f  %1.2f %1.2f\n",_7G80[44],_7G80[84],_7G80[125] );


    //5G
    _5G30(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5G30[i] = 2.82 + 0.00175*(i-5);
            else if (i<85 && i>=45) _5G30[i] = 2.89+0.0018*(i-45);
            else if (i<140 && i >=85) _5G30[i] = 2.96+0.0012*(i-85);
        }
        //printf("5G  %1.2f  %1.2f %1.2f\n",_5G30[44],_5G30[84],_5G30[125] );
    _5G40(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5G40[i] = 2.80 + 0.0015*(i-5);
            else if (i<85 && i>=45) _5G40[i] = 2.86+0.00175*(i-45);
            else if (i<140 && i >=85) _5G40[i] = 2.93+0.00125*(i-85);
        }
        //printf("5G  %1.2f  %1.2f %1.2f\n",_5G40[44],_5G40[84],_5G40[125] );
    _5G50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5G50[i] = 2.79 + 0.001*(i-5);
            else if (i<85 && i>=45) _5G50[i] = 2.84+0.0015*(i-45);
            else if (i<140 && i >=85) _5G50[i] = 2.90+0.0015*(i-85);
        }
        //printf("5G  %1.2f  %1.2f %1.2f\n",_5G50[44],_5G50[84],_5G50[125] );
    _5G60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5G60[i] = 2.78 + 0.001*(i-5);
            else if (i<85 && i>=45) _5G60[i] = 2.82+0.00175*(i-45);
            else if (i<140 && i >=85) _5G60[i] = 2.89+0.001*(i-85);
        }
        //printf("5G  %1.2f  %1.2f %1.2f\n",_5G60[44],_5G60[84],_5G60[125] );
    _5G70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5G70[i] = 2.77 + 0.001*(i-5);
            else if (i<85 && i>=45) _5G70[i] = 2.81+0.00125*(i-45);
            else if (i<140 && i >=85) _5G70[i] = 2.86+0.00125*(i-85);
        }
        //printf("5G  %1.2f  %1.2f %1.2f\n",_5G70[44],_5G70[84],_5G70[125] );
    _5G80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5G80[i] = 2.76 + 0.001*(i-5);
            else if (i<85 && i>=45) _5G80[i] = 2.8+0.00125*(i-45);
            else if (i<140 && i >=85) _5G80[i] = 2.85+0.00125*(i-85);
        }
        //printf("5G  %1.2f  %1.2f %1.2f\n",_5G80[44],_5G80[84],_5G80[125] );

    //25G
    _25G30(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _25G30[i] = 2.68 + 0.0015*(i-5);
            else if (i<85 && i>=45) _25G30[i] = 2.74+0.0018*(i-45);
            else if (i<140 && i >=85) _25G30[i] = 2.81+0.002*(i-85);
        }
        //printf("25G  %1.2f  %1.2f %1.2f\n",_25G30[44],_25G30[84],_25G30[125] );
    _25G40(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _25G40[i] = 2.68 + 0.00075*(i-5);
            else if (i<85 && i>=45) _25G40[i] = 2.71+0.0015*(i-45);
            else if (i<140 && i >=85) _25G40[i] = 2.77+0.00125*(i-85);
        }
        //printf("25G  %1.2f  %1.2f %1.2f\n",_25G40[44],_25G40[84],_25G40[125] );
    _25G50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _25G50[i] = 2.65 + 0.00075*(i-5);
            else if (i<85 && i>=45) _25G50[i] = 2.68+0.00125*(i-45);
            else if (i<140 && i >=85) _25G50[i] = 2.73+0.00125*(i-85);
        }
        //printf("25G  %1.2f  %1.2f %1.2f\n",_25G50[44],_25G50[84],_25G50[125] );
    _25G60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _25G60[i] = 2.64 + 0.0005*(i-5);
            else if (i<85 && i>=45) _25G60[i] = 2.66+0.001*(i-45);
            else if (i<140 && i >=85) _25G60[i] = 2.70+0.001*(i-85);
        }
        //printf("25G  %1.2f  %1.2f %1.2f\n",_25G60[44],_25G60[84],_25G60[125] );
    _25G70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _25G70[i] = 2.64 + 0.00*(i-5);
            else if (i<85 && i>=45) _25G70[i] = 2.64+0.00075*(i-45);
            else if (i<140 && i >=85) _25G70[i] = 2.67+0.001*(i-85);
        }
        //printf("25G  %1.2f  %1.2f %1.2f\n",_25G70[44],_25G70[84],_25G70[125] );
    _25G80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _25G80[i] = 2.63 + 0.00*(i-5);
            else if (i<85 && i>=45) _25G80[i] = 2.63+0.0005*(i-45);
            else if (i<140 && i >=85) _25G80[i] = 2.65+0.0005*(i-85);
        }
        //printf("25G  %1.2f  %1.2f %1.2f\n",_25G80[44],_25G80[84],_25G80[125] );


    //1G
    _1G30(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1G30[i] = 2.58 + 0.00025*(i-5);
            else if (i<85 && i>=45) _1G30[i] = 2.59+0.001*(i-45);
            else if (i<140 && i >=85) _1G30[i] = 2.63+0.00125*(i-85);
        }
        //printf("1G  %1.2f  %1.2f %1.2f\n",_1G30[44],_1G30[84],_1G30[125] );
    _1G40(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1G40[i] = 2.56 - 0.00025*(i-5);
            else if (i<85 && i>=45) _1G40[i] = 2.55+0.0005*(i-45);
            else if (i<140 && i >=85) _1G40[i] = 2.57+0.0005*(i-85);
        }
        //printf("1G  %1.2f  %1.2f %1.2f\n",_1G40[44],_1G40[84],_1G40[125] );
    _1G50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1G50[i] = 2.55 - 0.00025*(i-5);
            else if (i<85 && i>=45) _1G50[i] = 2.54+0.00025*(i-45);
            else if (i<140 && i >=85) _1G50[i] = 2.55+0.0005*(i-85);
        }
        //printf("1G  %1.2f  %1.2f %1.2f\n",_1G50[44],_1G50[84],_1G50[125] );
    _1G60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1G60[i] = 2.54 - 0.0005*(i-5);
            else if (i<85 && i>=45) _1G60[i] = 2.52+0.00025*(i-45);
            else if (i<140 && i >=85) _1G60[i] = 2.53+0.00025*(i-85);
        }
        //printf("1G  %1.2f  %1.2f %1.2f\n",_1G60[44],_1G60[84],_1G60[125] );
    _1G70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1G70[i] = 2.53 - 0.0005*(i-5);
            else if (i<85 && i>=45) _1G70[i] = 2.51+0.0*(i-45);
            else if (i<140 && i >=85) _1G70[i] = 2.51+0.00025*(i-85);
        }
        //printf("1G  %1.2f  %1.2f %1.2f\n",_1G70[44],_1G70[84],_1G70[125] );
    _1G80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _1G80[i] = 2.52 - 0.0005*(i-5);
            else if (i<85 && i>=45) _1G80[i] = 2.50+0.00*(i-45);
            else if (i<140 && i >=85) _1G80[i] = 2.50+0.00*(i-85);
        }
        //printf("1G  %1.2f  %1.2f %1.2f\n",_1G80[44],_1G80[84],_1G80[125] );


    //10GY
    _10GY30(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10GY30[i] = 2.52 - 0.001*(i-5);
            else if (i<85 && i>=45) _10GY30[i] = 2.48-0.002*(i-45);
            else if (i<140 && i >=85) _10GY30[i] = 2.40+0.0025*(i-85);
        }
        //printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY30[44],_10GY30[84],_10GY30[125] );
    _10GY40(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10GY40[i] = 2.48 - 0.0005*(i-5);
            else if (i<85 && i>=45) _10GY40[i] = 2.46-0.0005*(i-45);
            else if (i<140 && i >=85) _10GY40[i] = 2.44-0.0015*(i-85);
        }
        //printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY40[44],_10GY40[84],_10GY40[125] );
    _10GY50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10GY50[i] = 2.48 - 0.00075*(i-5);
            else if (i<85 && i>=45) _10GY50[i] = 2.45-0.00075*(i-45);
            else if (i<140 && i >=85) _10GY50[i] = 2.42-0.00175*(i-85);
        }
        //printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY50[44],_10GY50[84],_10GY50[125] );
    _10GY60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10GY60[i] = 2.47 - 0.00125*(i-5);
            else if (i<85 && i>=45) _10GY60[i] = 2.42-0.00025*(i-45);
            else if (i<140 && i >=85) _10GY60[i] = 2.41-0.0005*(i-85);
        }
        //printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY60[44],_10GY60[84],_10GY60[125] );
    _10GY70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10GY70[i] = 2.46 - 0.001*(i-5);
            else if (i<85 && i>=45) _10GY70[i] = 2.42+0.0*(i-45);
            else if (i<140 && i >=85) _10GY70[i] = 2.42-0.001*(i-85);
        }
        //printf("10GY %1.2f  %1.2f %1.2f\n",_10GY70[44],_10GY70[84],_10GY70[125] );
    _10GY80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _10GY80[i] = 2.45 - 0.00075*(i-5);
            else if (i<85 && i>=45) _10GY80[i] = 2.42 - 0.0005*(i-45);
            else if (i<140 && i >=85) _10GY80[i] = 2.40-0.0005*(i-85);
        }
        //printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY80[44],_10GY80[84],_10GY80[125] );


    //75GY
    _75GY30(maxInd2);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _75GY30[i] = 2.36 - 0.0025*(i-5);
            else if (i<90 && i>=45) _75GY30[i] = 2.26-0.00175*(i-45);
        }
        //printf("75GY  %1.2f  %1.2f\n",_75GY30[44],_75GY30[84] );
    _75GY40(maxInd2);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _75GY40[i] = 2.34 - 0.00175*(i-5);
            else if (i<90 && i>=45) _75GY40[i] = 2.27-0.00225*(i-45);
        }
        //printf("75GY  %1.2f  %1.2f \n",_75GY40[44],_75GY40[84] );
    _75GY50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _75GY50[i] = 2.32 - 0.0015*(i-5);
            else if (i<85 && i>=45) _75GY50[i] = 2.26-0.00175*(i-45);
            else if (i<140 && i >=85) _75GY50[i] = 2.19-0.00325*(i-85);
        }
        //printf("75GY  %1.2f  %1.2f %1.2f %1.2f\n",_75GY50[44],_75GY50[84],_75GY50[125],_75GY50[139] );
    _75GY60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _75GY60[i] = 2.30 - 0.00125*(i-5);
            else if (i<85 && i>=45) _75GY60[i] = 2.25-0.001*(i-45);
            else if (i<140 && i >=85) _75GY60[i] = 2.21-0.0027*(i-85);
        }
        //printf("75GY  %1.2f  %1.2f %1.2f\n",_75GY60[44],_75GY60[84],_75GY60[125] );
    _75GY70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _75GY70[i] = 2.29 - 0.00125*(i-5);
            else if (i<85 && i>=45) _75GY70[i] = 2.24-0.0015*(i-45);
            else if (i<140 && i >=85) _75GY70[i] = 2.18-0.00175*(i-85);
        }
        //printf("75GY %1.2f  %1.2f %1.2f\n",_75GY70[44],_75GY70[84],_75GY70[125] );
    _75GY80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _75GY80[i] = 2.27 - 0.001*(i-5);
            else if (i<85 && i>=45) _75GY80[i] = 2.23 - 0.001*(i-45);
            else if (i<140 && i >=85) _75GY80[i] = 2.19-0.00175*(i-85);
        }
        //printf("75GY  %1.2f  %1.2f %1.2f\n",_75GY80[44],_75GY80[84],_75GY80[125] );


    //55GY
    _5GY30(maxInd2);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5GY30[i] = 2.16 - 0.002*(i-5);
            else if (i<90 && i>=45) _5GY30[i] = 2.07-0.0025*(i-45);
        }
        //printf("5GY  %1.2f  %1.2f\n",_5GY30[44],_5GY30[84] );

    //5GY4: 2.14,2.04, 1.96, 1.91 //95

    _5GY40(maxInd2);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5GY40[i] = 2.14 - 0.0025*(i-5);
            else if (i<90 && i>=45) _5GY40[i] = 2.04-0.003*(i-45);
        }
        //printf("5GY  %1.2f  %1.2f \n",_5GY40[44],_5GY40[84] );
    _5GY50(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5GY50[i] = 2.13 - 0.00175*(i-5);
            else if (i<85 && i>=45) _5GY50[i] = 2.06-0.002*(i-45);
            else if (i<140 && i >=85) _5GY50[i] = 1.98-0.00225*(i-85);
        }
        //printf("5GY  %1.2f  %1.2f %1.2f\n",_5GY50[44],_5GY50[84],_5GY50[125] );
    _5GY60(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5GY60[i] = 2.11 - 0.0015*(i-5);
            else if (i<85 && i>=45) _5GY60[i] = 2.05-0.002*(i-45);
            else if (i<140 && i >=85) _5GY60[i] = 1.97-0.00275*(i-85);
        }
        //printf("5GY  %1.2f  %1.2f %1.2f\n",_5GY60[44],_5GY60[84],_5GY60[125] );
    _5GY70(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5GY70[i] = 2.09 - 0.001*(i-5);
            else if (i<85 && i>=45) _5GY70[i] = 2.05-0.00175*(i-45);
            else if (i<140 && i >=85) _5GY70[i] = 1.98-0.002*(i-85);
        }
        //printf("5GY %1.2f  %1.2f %1.2f\n",_5GY70[44],_5GY70[84],_5GY70[125] );
    _5GY80(maxInd);
        for (int i=0; i<maxInd; i++) {
            if (i<45 && i>5) _5GY80[i] = 2.07 - 0.001*(i-5);
            else if (i<85 && i>=45) _5GY80[i] = 2.03 - 0.00075*(i-45);
            else if (i<140 && i >=85) _5GY80[i] = 2.0-0.002*(i-85);
        }
        //printf("5GY  %1.2f  %1.2f %1.2f\n",_5GY80[44],_5GY80[84],_5GY80[125] );

    #ifdef _DEBUG
        t2e.set();
        if (settings->verbose)
            printf("Lutf Munsell  %d usec\n", t2e.etime(t1e));
    #endif
    }

}
