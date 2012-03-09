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
#include <cmath>
#include "../rtgui/mycurve.h"
#include "../rtgui/myflatcurve.h"
#include "../rtgui/mydiagonalcurve.h"

#include "LUT.h"

#define CURVES_MIN_POLY_POINTS  1000

#define SQR(x) ((x)*(x))

#define CLIPI(a) ((a)>0?((a)<65534?(a):65534):0)

namespace rtengine {

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
        double lx = log(x);
        return m2*x + (1.0-m2)*(2.0 - exp(k*lx))*exp(l*lx);
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

    static void init ();
    static void cleanup ();

    static inline double centercontrast   (double x, double b, double m);
    
	static inline float hlcurve (const float exp_scale, const float comp, const float hlrange, float level) 
	{
		if (comp>0.0) {
			float val = level+(hlrange-65536.0);
			float Y = val*exp_scale/hlrange;
			float R = hlrange/(val*comp);
			return log(1.0+Y*comp)*R;
		} else {
			return exp_scale;
		}
	}

  public:
    static void complexCurve (double ecomp, double black, double hlcompr, double hlcomprthresh, double shcompr, double br, double contr,
							  double gamma_, bool igamma_, const std::vector<double>& curvePoints, LUTu & histogram, LUTu & histogramCropped,
							  LUTf & hlCurve, LUTf & shCurve,LUTf & outCurve, LUTu & outBeforeCCurveHistogram, int skip=1);
	static void complexsgnCurve (double saturation, bool satlimit, double satlimthresh, const std::vector<double>& acurvePoints,
								 const std::vector<double>& bcurvePoints, LUTf & aoutCurve, LUTf & boutCurve, LUTf & satCurve, int skip=1);
	static void complexLCurve (double br, double contr, const std::vector<double>& curvePoints, LUTu & histogram, LUTu & histogramCropped,
							   LUTf & outCurve, LUTu & outBeforeCCurveHistogram, int skip); 
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
    void AddPolygons ();
    virtual double getVal (double t) = 0;
    virtual void   getVal (const std::vector<double>& t, std::vector<double>& res) = 0;

    virtual bool   isIdentity () = 0;
};

class DiagonalCurve : public Curve {

  protected:
    DiagonalCurveType kind;

    unsigned int minSearch;     // a effacer!!!
    unsigned int maxSearch;     // a effacer!!!
    unsigned int searchArray[21];     // a effacer!!!

    void spline_cubic_set ();
    void NURBS_set ();

  public:
    DiagonalCurve (const std::vector<double>& points, int ppn=CURVES_MIN_POLY_POINTS);
    virtual ~DiagonalCurve ();

    double getVal     (double t);
    void   getVal     (const std::vector<double>& t, std::vector<double>& res);
    bool   isIdentity () { return kind==DCT_Empty; };
};

class FlatCurve : public Curve {

  protected:
    FlatCurveType kind;
    double* leftTangent;
    double* rightTangent;
    bool periodic;

    void CtrlPoints_set ();

  public:

    FlatCurve (const std::vector<double>& points, bool isPeriodic = true, int ppn=CURVES_MIN_POLY_POINTS);
    virtual ~FlatCurve ();

    double getVal     (double t);
    void   getVal     (const std::vector<double>& t, std::vector<double>& res);
    bool   isIdentity () { return kind==FCT_Empty; };
};
}

#undef CLIPI

#endif
