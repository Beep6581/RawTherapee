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
#ifndef __CURVES_H__
#define __CURVES_H__

#include <glibmm.h>
#include <map>
#include <string>
#include "rt_math.h"
#include "../rtgui/mycurve.h"
#include "../rtgui/myflatcurve.h"
#include "../rtgui/mydiagonalcurve.h"
#include "color.h"
#include "procparams.h"
#include "editbuffer.h"

#include "LUT.h"

#define CURVES_MIN_POLY_POINTS  1000

#include "rt_math.h"

#define CLIPI(a) ((a)>0?((a)<65534?(a):65534):0)

using namespace std;

namespace rtengine {
    class ToneCurve;
    class ColorAppearance;

class CurveFactory {

    friend class Curve;

  protected:

    // functions calculating the parameters of the contrast curve based on the desired slope at the center
    static double solve_upper (double m, double c, double deriv);
    static double solve_lower (double m, double c, double deriv);
    static double dupper (const double b, const double m, const double c);
    static double dlower (const double b, const double m, const double c);
    
    // basic convex function between (0,0) and (1,1). m1 and m2 controls the slope at the start and end point
    static inline double basel (double x, double m1, double m2) {
        if (x==0.0)
            return 0.0;
        double k = sqrt ((m1-1.0)*(m1-m2)/2) / (1.0-m2);
        double l = (m1-m2) / (1.0-m2) + k;
        double lx = xlog(x);
        return m2*x + (1.0-m2)*(2.0 - xexp(k*lx))*xexp(l*lx);
    }
    // basic concave function between (0,0) and (1,1). m1 and m2 controls the slope at the start and end point
    static inline double baseu (double x, double m1, double m2) {
        return 1.0 - basel(1.0-x, m1, m2);
    }
    // convex curve between (0,0) and (1,1) with slope m at (0,0). hr controls the highlight recovery 
    static inline double cupper (double x, double m, double hr) {
        if (hr>1.0)
            return baseu (x, m, 2.0*(hr-1.0)/m);
        double x1 = (1.0-hr)/m;
        double x2 = x1 + hr;
        if (x>=x2) return 1.0;
        if (x<x1) return x*m;
        return 1.0 - hr + hr*baseu((x-x1)/hr, m, 0);
    }
    // concave curve between (0,0) and (1,1) with slope m at (1,1). sr controls the shadow recovery 
    static inline double clower (double x, double m, double sr) {
        return 1.0 - cupper(1.0-x, m, sr);
    }
	// convex curve between (0,0) and (1,1) with slope m at (0,0). hr controls the highlight recovery 
    static inline double cupper2 (double x, double m, double hr) {
        double x1 = (1.0-hr)/m;
        double x2 = x1 + hr;
        if (x>=x2) return 1.0;
        if (x<x1) return x*m;
        return 1.0 - hr + hr*baseu((x-x1)/hr, m, 0.3*hr);
    }
	static inline double clower2 (double x, double m, double sr) {
		//curve for b<0; starts with positive slope and then rolls over toward straight line to x=y=1
		float x1 = sr/1.5 + 0.00001;
		float y1 = 1-(1-x1)*m;
		if (x>x1 || sr<0.001) 
			return 1-(1-x)*m;
		else
			return y1+m*(x-x1)-(1-m)*SQR(SQR(1-x/x1));
	}
    // tone curve base. a: slope (from exp.comp.), b: black point normalized by 65535, 
	// D: max. x value (can be>1), hr,sr: highlight,shadow recovery
    static inline double basecurve (double x, double a, double b, double D, double hr, double sr) { 
        if (b<0) {
			double m = 0.5;//midpoint
			double slope = 1.0+b;//slope of straight line between (0,-b) and (1,1)
			double y = -b+m*slope;//value at midpoint
			if (x>m) 
				return y + (x - m)*slope;//value on straight line between (m,y) and (1,1)
			else 
				return y*clower2(x/m, slope*m/y, 2.0-sr);
		} else {
			double slope = a/(1.0-b);
			double m = a*D>1.0 ? b/a+(0.25)/slope : b+(1-b)/4;
			double y = a*D>1.0 ? 0.25 : (m-b/a)*slope;
			if (x<=m)
				return b==0 ? x*slope : clower (x/m, slope*m/y, sr) * y;
			else if (a*D>1.0)
				return y+(1.0-y)*cupper2((x-m)/(D-m), slope*(D-m)/(1.0-y), hr);
			else
				return y+(x-m)*slope;
		}
    }
	

  public:
	const static double sRGBGamma;  // standard average gamma
    const static double sRGBGammaCurve;  // 2.4 in the curve

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// accurately determine value from integer array with float as index
	//linearly interpolate from ends of range if arg is out of bounds
	static inline float interp(int *array,float f)
	{
		int index = CLIPI(floor(f));
		float part = (float)((f)-index)*(float)(array[index+1]-array[index]);
		return (float)array[index]+part;
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// accurately determine value from float array with float as index
	//linearly interpolate from ends of range if arg is out of bounds
	static inline float flinterp(float *array,float f)
	{
		int index = CLIPI(floor(f));
		float part = ((f)-(float)index)*(array[index+1]-array[index]);
		return array[index]+part;
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    static inline double centercontrast   (double x, double b, double m);
    
    // standard srgb gamma and its inverse
    static inline double gamma2            (double x) {
                                            return x <= 0.00304 ? x*12.92 : 1.055*exp(log(x)/sRGBGammaCurve)-0.055;
                                          }
    static inline double igamma2           (double x) {
                                            return x <= 0.03928 ? x/12.92 : exp(log((x+0.055)/1.055)*sRGBGammaCurve);
                                          }
    // gamma function with adjustable parameters
    static inline double gamma            (double x, double gamma, double start, double slope, double mul, double add){
                                            return (x <= start ? x*slope : exp(log(x)/gamma)*mul-add);
                                          }
    static inline double igamma           (double x, double gamma, double start, double slope, double mul, double add){
											return (x <= start*slope ? x/slope : exp(log((x+add)/mul)*gamma) );
										  }
											
	static inline float hlcurve (const float exp_scale, const float comp, const float hlrange, float level) 
	{
		if (comp>0.0) {
			float val = level+(hlrange-65536.0);
			if(val == 0.0f)	// to avoid division by zero
				val = 0.000001f;
			float Y = val*exp_scale/hlrange;
			Y *= comp;
			if(Y <= -1.0)	// to avoid log(<=0)
				Y = -.999999f;
			float R = hlrange/(val*comp);
			return log(1.0+Y)*R;
		} else {
			return exp_scale;
		}
	}

  public:
    static void complexCurve (double ecomp, double black, double hlcompr, double hlcomprthresh, double shcompr, double br, double contr,
							  double gamma_, bool igamma_, procparams::ToneCurveParams::eTCModeId curveMode, const std::vector<double>& curvePoints, procparams::ToneCurveParams::eTCModeId curveMode2, const std::vector<double>& curvePoints2, 
							  
							  LUTu & histogram, LUTu & histogramCropped,
							  LUTf & hlCurve, LUTf & shCurve,LUTf & outCurve, LUTu & outBeforeCCurveHistogram, ToneCurve & outToneCurve, ToneCurve & outToneCurve2, 
				  
							  int skip=1);
	static void curveBW (const std::vector<double>& curvePointsbw, const std::vector<double>& curvePointsbw2, LUTu & histogrambw, LUTu & outBeforeCCurveHistogrambw,
						 ToneCurve & customToneCurvebw1, ToneCurve & customToneCurvebw2, int skip);
		
	static void curveCL ( bool & clcutili, const std::vector<double>& clcurvePoints, LUTf & clCurve, LUTu & histogramcl, LUTu & outBeforeCLurveHistogram, int skip);
	static void curveToningCL ( bool & clctoningutili, const std::vector<double>& clcurvePoints, LUTf & clToningCurve, int skip);
	static void curveToningLL ( bool & llctoningutili, const std::vector<double>& llcurvePoints, LUTf & llToningCurve, int skip);
						  
	static void complexsgnCurve ( float adjustr, bool & autili,  bool & butili, bool & ccutili, bool & clcutili, double saturation, double rstprotection, const std::vector<double>& acurvePoints,
								 const std::vector<double>& bcurvePoints,const std::vector<double>& cccurvePoints,const std::vector<double>& lccurvePoints, LUTf & aoutCurve, LUTf & boutCurve, LUTf & satCurve, LUTf & lhskCurve, 
								 LUTu & histogramC, LUTu & histogramLC, LUTu & outBeforeCCurveHistogram,LUTu & outBeforeLCurveHistogram,///for chroma
								int skip=1);
	static void complexLCurve (double br, double contr, const std::vector<double>& curvePoints, LUTu & histogram, LUTu & histogramCropped,
							   LUTf & outCurve, LUTu & outBeforeCCurveHistogram, int skip, bool & utili); 

	static void updatechroma (
					const std::vector<double>& cccurvePoints,
					LUTu & histogramC, LUTu & outBeforeCCurveHistogramC,//for chroma
					int skip=1);

	static void curveLightBrightColor (
					procparams::ColorAppearanceParams::eTCModeId curveMode, const std::vector<double>& curvePoints,
					procparams::ColorAppearanceParams::eTCModeId curveMode2, const std::vector<double>& curvePoints2,
					procparams::ColorAppearanceParams::eCTCModeId curveMode3, const std::vector<double>& curvePoints3,
					LUTu & histogram, LUTu & histogramCropped, LUTu & outBeforeCCurveHistogram,
					LUTu & histogramC, LUTu & outBeforeCCurveHistogramC,
					ColorAppearance & outColCurve1,
					ColorAppearance & outColCurve2,
					ColorAppearance & outColCurve3,
					int skip=1);
	static void RGBCurve (const std::vector<double>& curvePoints, LUTf & outCurve, int skip);

};

class Curve {

    class HashEntry {
    public:
        unsigned short smallerValue;
        unsigned short higherValue;
    };
  protected:
    int N;
    int ppn;			// targeted polyline point number
    double* x;
    double* y;
    // begin of variables used in Parametric curves only
    double mc;
    double mfc;
    double msc;
    double mhc;
    // end of variables used in Parametric curves only
    std::vector<double> poly_x;		// X points of the faceted curve
    std::vector<double> poly_y;		// Y points of the faceted curve
    std::vector<HashEntry> hash;
    unsigned short hashSize;		// hash table's size, between [10, 100, 1000]

    double* ypp;

    // Fields for the elementary curve polygonisation
    double x1, y1, x2, y2, x3, y3;
    bool firstPointIncluded;
    double increment;
    int nbr_points;

    static inline double p00 (double x, double prot) { return CurveFactory::clower (x, 2.0, prot); }
    static inline double p11 (double x, double prot) { return CurveFactory::cupper (x, 2.0, prot); }
    static inline double p01 (double x, double prot) { return x<=0.5 ? CurveFactory::clower (x*2, 2.0, prot)/2.0 : 0.5 + CurveFactory::cupper ((x-0.5)*2, 2.0, prot)/2.0; }
    static inline double p10 (double x, double prot) { return x<=0.5 ? CurveFactory::cupper (x*2, 2.0, prot)/2.0 : 0.5 + CurveFactory::clower ((x-0.5)*2, 2.0, prot)/2.0; }
    static inline double pfull (double x, double prot, double sh, double hl) { return (1-sh)*(1-hl)*p00(x,prot) + sh*hl*p11(x,prot) + (1-sh)*hl*p01(x,prot) + sh*(1-hl)*p10(x,prot); }

    void fillHash();

  public:
    Curve ();
    virtual ~Curve () {};
    void AddPolygons ();
    int getSize () const; // return the number of control points
    void getControlPoint(int cpNum, double &x, double &y) const;
    virtual double getVal (double t) const = 0;
    virtual void   getVal (const std::vector<double>& t, std::vector<double>& res) const = 0;

    virtual bool   isIdentity () const = 0;
};

class DiagonalCurve : public Curve {

  protected:
    DiagonalCurveType kind;

    void spline_cubic_set ();
    void NURBS_set ();

  public:
    DiagonalCurve (const std::vector<double>& points, int ppn=CURVES_MIN_POLY_POINTS);
    virtual ~DiagonalCurve ();

    double getVal     (double t) const;
    void   getVal     (const std::vector<double>& t, std::vector<double>& res) const;
    bool   isIdentity () const { return kind==DCT_Empty; };
};

class FlatCurve : public Curve {

  protected:
    FlatCurveType kind;
    double* leftTangent;
    double* rightTangent;
    double identityValue;
    bool periodic;

    void CtrlPoints_set ();

  public:

    FlatCurve (const std::vector<double>& points, bool isPeriodic = true, int ppn=CURVES_MIN_POLY_POINTS);
    virtual ~FlatCurve ();

    double getVal     (double t) const;
    void   getVal     (const std::vector<double>& t, std::vector<double>& res) const;
    bool   setIdentityValue (double iVal);
    bool   isIdentity () const { return kind==FCT_Empty; };
};

class ToneCurve {
  public:
    LUTf lutToneCurve;  // 0xffff range

    virtual ~ToneCurve() {};

    void Reset();
    void Set(Curve *pCurve);
    operator bool (void) const { return lutToneCurve; }
};

class OpacityCurve {
  public:
    LUTf lutOpacityCurve;  // 0xffff range

    virtual ~OpacityCurve() {};

    void Reset();
    void Set(const Curve *pCurve);
    void Set(const std::vector<double> &curvePoints);

    // TODO: transfer this method to the Color class...
    float blend (float x, float lower, float upper) const {
        return (upper-lower)*lutOpacityCurve[x*500.f] + lower;
    }
    void blend3f (float x, float lower1, float upper1, float &result1, float lower2, float upper2, float &result2, float lower3, float upper3, float &result3) const {
        float opacity = lutOpacityCurve[x*500.f];
        result1 = (upper1-lower1)*opacity + lower1;
        result2 = (upper2-lower2)*opacity + lower2;
        result3 = (upper3-lower3)*opacity + lower3;
    }

    operator bool (void) const { return lutOpacityCurve; }
};

class ColorGradientCurve {
  public:
    LUTf   lut1;    // [0.;1.] range (float values)
    LUTf   lut2;  // [0.;1.] range (float values)
    LUTf   lut3;   // [0.;1.] range (float values)
    double low;
    double high;

    virtual ~ColorGradientCurve() {};

    void Reset();
    void SetXYZ(const Curve *pCurve, const double xyz_rgb[3][3], const double rgb_xyz[3][3], float satur, float lumin);
    void SetXYZ(const std::vector<double> &curvePoints, const double xyz_rgb[3][3], const double rgb_xyz[3][3], float satur, float lumin);
    void SetRGB(const Curve *pCurve, const double xyz_rgb[3][3], const double rgb_xyz[3][3]);
    void SetRGB(const std::vector<double> &curvePoints, const double xyz_rgb[3][3], const double rgb_xyz[3][3]);

    /**
    * @brief Get the value of Red, Green and Blue corresponding to the requested index
    * @param index value in the [0 ; 1] range
    * @param r corresponding red value [0 ; 65535] (return value)
    * @param g corresponding green value [0 ; 65535] (return value)
    * @param b corresponding blue value [0 ; 65535] (return value)
    */
    void getVal(float index, float &r, float &g, float &b) const;
    operator bool (void) const { return lut1 && lut2 && lut3; }
};

class ColorAppearance {
  public:
    LUTf lutColCurve;  // 0xffff range

    virtual ~ColorAppearance() {};

    void Reset();
    void Set(Curve *pCurve);
    operator bool (void) const { return lutColCurve; }
};

class Lightcurve : public ColorAppearance {
  public:
    void Apply(float& Li) const;
};

//lightness curve
inline void Lightcurve::Apply (float& Li) const {

    assert (lutColCurve);

    Li = lutColCurve[Li];
}

class Brightcurve : public ColorAppearance {
  public:
    void Apply(float& Br) const;
};

//brightness curve
inline void Brightcurve::Apply (float& Br) const {

    assert (lutColCurve);

    Br = lutColCurve[Br];
}

class Chromacurve : public ColorAppearance {
  public:
    void Apply(float& Cr) const;
};

//Chroma curve
inline void Chromacurve::Apply (float& Cr) const {

    assert (lutColCurve);

    Cr = lutColCurve[Cr];
}
class Saturcurve : public ColorAppearance {
  public:
    void Apply(float& Sa) const;
};

//Saturation curve
inline void Saturcurve::Apply (float& Sa) const {

    assert (lutColCurve);

    Sa = lutColCurve[Sa];
}

class Colorfcurve : public ColorAppearance {
  public:
    void Apply(float& Cf) const;
};

//Colorfullness curve
inline void Colorfcurve::Apply (float& Cf) const {

    assert (lutColCurve);

    Cf = lutColCurve[Cf];
}


class StandardToneCurve : public ToneCurve {
  public:
    void Apply(float& r, float& g, float& b) const;
};
class StandardToneCurvebw : public ToneCurve {
  public:
    void Apply(float& r, float& g, float& b) const;
};

class AdobeToneCurve : public ToneCurve {
  private:
    void RGBTone(float& r, float& g, float& b) const;  // helper for tone curve

  public:
    void Apply(float& r, float& g, float& b) const;
};

class AdobeToneCurvebw : public ToneCurve {
  private:
    void RGBTone(float& r, float& g, float& b) const;  // helper for tone curve

  public:
    void Apply(float& r, float& g, float& b) const;
};

class SatAndValueBlendingToneCurve : public ToneCurve {
  public:
    void Apply(float& r, float& g, float& b) const;
};

class SatAndValueBlendingToneCurvebw : public ToneCurve {
  public:
    void Apply(float& r, float& g, float& b) const;
};

class WeightedStdToneCurve : public ToneCurve {
  private:
    float Triangle(float refX, float refY, float X2) const;
  public:
    void Apply(float& r, float& g, float& b) const;
};

class WeightedStdToneCurvebw : public ToneCurve {
  private:
    float Triangle(float refX, float refY, float X2) const;
  public:
    void Apply(float& r, float& g, float& b) const;
};

// Standard tone curve
inline void StandardToneCurve::Apply (float& r, float& g, float& b) const {

    assert (lutToneCurve);

    r = lutToneCurve[r];
    g = lutToneCurve[g];
    b = lutToneCurve[b];
}
// Standard tone curve
inline void StandardToneCurvebw::Apply (float& r, float& g, float& b) const {

    assert (lutToneCurve);

    r = lutToneCurve[r];
    g = lutToneCurve[g];
    b = lutToneCurve[b];
}

// Tone curve according to Adobe's reference implementation
// values in 0xffff space
// inlined to make sure there will be no cache flush when used
inline void AdobeToneCurve::Apply (float& r, float& g, float& b) const {

    assert (lutToneCurve);

    if (r >= g) {
        if      (g > b) RGBTone (r, g, b); // Case 1: r >= g >  b
        else if (b > r) RGBTone (b, r, g); // Case 2: b >  r >= g
        else if (b > g) RGBTone (r, b, g); // Case 3: r >= b >  g
        else {                             // Case 4: r >= g == b
            r = lutToneCurve[r];
            g = lutToneCurve[g];
            b = g;
        }
    }
    else {
        if      (r >= b) RGBTone (g, r, b); // Case 5: g >  r >= b
        else if (b >  g) RGBTone (b, g, r); // Case 6: b >  g >  r
        else             RGBTone (g, b, r); // Case 7: g >= b >  r
    }
}
inline void AdobeToneCurvebw::Apply (float& r, float& g, float& b) const {

    assert (lutToneCurve);

    if (r >= g) {
        if      (g > b) RGBTone (r, g, b); // Case 1: r >= g >  b
        else if (b > r) RGBTone (b, r, g); // Case 2: b >  r >= g
        else if (b > g) RGBTone (r, b, g); // Case 3: r >= b >  g
        else {                             // Case 4: r >= g == b
            r = lutToneCurve[r];
            g = lutToneCurve[g];
            b = g;
        }
    }
    else {
        if      (r >= b) RGBTone (g, r, b); // Case 5: g >  r >= b
        else if (b >  g) RGBTone (b, g, r); // Case 6: b >  g >  r
        else             RGBTone (g, b, r); // Case 7: g >= b >  r
    }
}

inline void AdobeToneCurve::RGBTone (float& r, float& g, float& b) const {
    float rold=r,gold=g,bold=b;

    r = lutToneCurve[rold];
    b = lutToneCurve[bold];
    g = b + ((r - b) * (gold - bold) / (rold - bold));
}
inline void AdobeToneCurvebw::RGBTone (float& r, float& g, float& b) const {
    float rold=r,gold=g,bold=b;

    r = lutToneCurve[rold];
    b = lutToneCurve[bold];
    g = b + ((r - b) * (gold - bold) / (rold - bold));
}

inline float WeightedStdToneCurve::Triangle(float a, float a1, float b) const {
	if (a != b) {
		float b1;
		float a2 = a1 - a;
		if (b < a) { b1 = b + a2 *      b  /     a ; }
		else       { b1 = b + a2 * (65535.f-b) / (65535.f-a); }
		return b1;
	}
	return a1;
}
inline float WeightedStdToneCurvebw::Triangle(float a, float a1, float b) const {
	if (a != b) {
		float b1;
		float a2 = a1 - a;
		if (b < a) { b1 = b + a2 *      b  /     a ; }
		else       { b1 = b + a2 * (65535.f-b) / (65535.f-a); }
		return b1;
	}
	return a1;
}

// Tone curve modifying the value channel only, preserving hue and saturation
// values in 0xffff space
inline void WeightedStdToneCurve::Apply (float& r, float& g, float& b) const {

    assert (lutToneCurve);

    float r1 = lutToneCurve[r];
    float g1 = Triangle(r, r1, g);
    float b1 = Triangle(r, r1, b);

    float g2 = lutToneCurve[g];
    float r2 = Triangle(g, g2, r);
    float b2 = Triangle(g, g2, b);

    float b3 = lutToneCurve[b];
    float r3 = Triangle(b, b3, r);
    float g3 = Triangle(b, b3, g);

    r = CLIP<float>( r1*0.50f + r2*0.25f + r3*0.25f);
    g = CLIP<float>(g1*0.25f + g2*0.50f + g3*0.25f);
    b = CLIP<float>(b1*0.25f + b2*0.25f + b3*0.50f);
}

inline void WeightedStdToneCurvebw::Apply (float& r, float& g, float& b) const {

    assert (lutToneCurve);

    float r1 = lutToneCurve[r];
    float g1 = Triangle(r, r1, g);
    float b1 = Triangle(r, r1, b);

    float g2 = lutToneCurve[g];
    float r2 = Triangle(g, g2, r);
    float b2 = Triangle(g, g2, b);

    float b3 = lutToneCurve[b];
    float r3 = Triangle(b, b3, r);
    float g3 = Triangle(b, b3, g);

    r = CLIP<float>( r1*0.50f + r2*0.25f + r3*0.25f);
    g = CLIP<float>(g1*0.25f + g2*0.50f + g3*0.25f);
    b = CLIP<float>(b1*0.25f + b2*0.25f + b3*0.50f);
}

// Tone curve modifying the value channel only, preserving hue and saturation
// values in 0xffff space
inline void SatAndValueBlendingToneCurve::Apply (float& r, float& g, float& b) const {

    assert (lutToneCurve);

    float h, s, v;
    float lum = (r+g+b)/3.f;
    //float lum = Color::rgbLuminance(r, g, b);
    float newLum = lutToneCurve[lum];
    if (newLum == lum)
    	return;
    bool increase = newLum > lum;

    Color::rgb2hsv(r, g, b, h, s, v);

    if (increase) {
    	// Linearly targeting Value = 1 and Saturation = 0
        float coef = (newLum-lum)/(65535.f-lum);
        float dV = (1.f-v)*coef;
        s *= 1.f-coef;
    	Color::hsv2rgb(h, s, v+dV, r, g, b);
    }
    else {
    	// Linearly targeting Value = 0
        float coef = (lum-newLum)/lum ;
        float dV = v*coef;
    	Color::hsv2rgb(h, s, v-dV, r, g, b);
    }
}

inline void SatAndValueBlendingToneCurvebw::Apply (float& r, float& g, float& b) const {

    assert (lutToneCurve);

    float h, s, v;
    float lum = (r+g+b)/3.f;
    //float lum = Color::rgbLuminance(r, g, b);
    float newLum = lutToneCurve[lum];
    if (newLum == lum)
    	return;
    bool increase = newLum > lum;

    Color::rgb2hsv(r, g, b, h, s, v);

    if (increase) {
    	// Linearly targeting Value = 1 and Saturation = 0
        float coef = (newLum-lum)/(65535.f-lum);
        float dV = (1.f-v)*coef;
        s *= 1.f-coef;
    	Color::hsv2rgb(h, s, v+dV, r, g, b);
    }
    else {
    	// Linearly targeting Value = 0
        float coef = (lum-newLum)/lum ;
        float dV = v*coef;
    	Color::hsv2rgb(h, s, v-dV, r, g, b);
    }
}

}

#undef CLIPI

#endif
