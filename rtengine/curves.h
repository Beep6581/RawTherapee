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

#include <math.h>
#include "rtcommon.h"

namespace rtengine {

class CurveFactory {

    friend class Curve;

  protected:

    // look-up tables for the standard srgb gamma and its inverse (filled by init())
    static int igammatab_srgb[65536];
    static int gammatab_srgb[65536];
    // look-up tables for the simple exponential gamma
    static int gammatab[65536];
    
    // functions calculating the parameters of the contrast curve based on the desired slope at the center
    static float solve_upper (float m, float c, float deriv);
    static float solve_lower (float m, float c, float deriv);
    static float dupper (const float b, const float m, const float c);
    static float dlower (const float b, const float m, const float c);
    
    // basic convex function between (0,0) and (1,1). m1 and m2 controls the slope at the start and end point
    static inline float basel (float x, float m1, float m2) {
        if (x==0.0)
            return 0.0;
        float k = sqrt ((m1-1.0)*(m1-m2)/2) / (1.0-m2);
        float l = (m1-m2) / (1.0-m2) + k;
        float lx = log(x);
        return m2*x + (1.0-m2)*(2.0 - exp(k*lx))*exp(l*lx);
    }
    // basic concave function between (0,0) and (1,1). m1 and m2 controls the slope at the start and end point
    static inline float baseu (float x, float m1, float m2) {
        return 1.0 - basel(1.0-x, m1, m2);
    }
    // convex curve between (0,0) and (1,1) with slope m at (0,0). hr controls the highlight recovery 
    static inline float cupper (float x, float m, float hr) {
        if (hr>1.0)
            return baseu (x, m, 2.0*(hr-1.0)/m);
        float x1 = (1.0-hr)/m;
        float x2 = x1 + hr;
        if (x>=x2) return 1.0;
        if (x<x1) return x*m;
        return 1.0 - hr + hr*baseu((x-x1)/hr, m, 0);
    }
    // concave curve between (0,0) and (1,1) with slope m at (1,1). sr controls the shadow recovery 
    static inline float clower (float x, float m, float sr) {
        return 1.0 - cupper(1.0-x, m, sr);
    }
    // tone curve base. a: slope (from exp.comp.), b: black, D: max. x value (can be>1), hr,sr: highlight,shadow recovery
    static inline float basecurve (float x, float a, float b, float D, float hr, float sr) {
        float m = b+0.5/a<D ? b+0.5/a : D;
        float y = (D-b)*a<0.5 ? (D-b)*a : 0.5;
        if (x<=m)
            return b==0 ? x*a : clower (x/m, a*m/y, sr) * y;
        else if (b+1.0/a<D)
            return y+(1.0-y)*cupper((x-m)/(D-m), a*(D-m)/(1.0-y), hr);
        else
            return y+(x-m)*a;
    }
    // brightness curve at point x, only positive amount it supported
    static inline float brightnessbase (float x, float amount) {
        if (x<0.5)
            return x + amount*cupper(2.0*x, 4.5, 0.0)/3.0;
        else
            return x + amount*cupper(2.0-2.0*x, 1.5, 0.0)/3.0;
    }
    // brightness curve at point x, positive negative and zero amount are supported
    static inline float brightness (float x, float amount) {
        if (amount==0)
            return x;
        else if (amount>0)
            return brightnessbase (x, amount);
        else 
            return 1.0 - brightnessbase (1.0-x, -amount);
    }

  public:

    static void init ();

    static inline float centercontrast   (float x, float b, float m);
    
    // standard srgb gamma and its inverse
    static inline float gamma2            (float x) {
                                            return x <= 0.00304 ? x*12.92 : 1.055*exp(log(x)/2.4)-0.055;
                                          }
    static inline float igamma2           (float x) {
                                            return x <= 0.03928 ? x/12.92 : exp(log((x+0.055)/1.055)*2.4);
                                          }
    // gamma function with adjustable parameters
    static inline float gamma            (float x, float gamma, float start, float slope, float mul, float add){
                                            return (x <= start ? x*slope : exp(log(x)/gamma)*mul-add);
                                          }

    // gamma functions on [0,65535] based on look-up tables
    static inline int    gamma_srgb       (int x) { return gammatab_srgb[x]; }
    static inline int    gamma            (int x) { return gammatab[x]; }
    static inline int    igamma_srgb      (int x) { return igammatab_srgb[x]; }

  public:
//    static void updateCurve3 (int* curve, int* ohistogram, const std::vector<float>& cpoints, float defmul, float ecomp, int black, float hlcompr, float shcompr, float br, float contr, float gamma_, bool igamma, int skip=1);
    static void complexCurve (float ecomp, float black, float hlcompr, float shcompr, float br, float contr, float gamma_, bool igamma, const FloatList& curvePoints, unsigned int* histogram, float* outCurve, unsigned int* outBeforeCCurveHistogram, int skip=1);

};

class Curve {

  protected:
    int N;
    float* x;
    float* y;
    float* ypp;
    int kind;       // = -1: linear interp., 0: empty, 1: spline interp., 2: parametric

  protected:
    void spline_cubic_set ();
    static inline float p00 (float x, float prot) { return CurveFactory::clower (x, 2.0, prot); }
    static inline float p11 (float x, float prot) { return CurveFactory::cupper (x, 2.0, prot); }
    static inline float p01 (float x, float prot) { return x<=0.5 ? CurveFactory::clower (x*2, 2.0, prot)/2.0 : 0.5 + CurveFactory::cupper ((x-0.5)*2, 2.0, prot)/2.0; }
    static inline float p10 (float x, float prot) { return x<=0.5 ? CurveFactory::cupper (x*2, 2.0, prot)/2.0 : 0.5 + CurveFactory::clower ((x-0.5)*2, 2.0, prot)/2.0; }
    static inline float pfull (float x, float prot, float sh, float hl) { return (1-sh)*(1-hl)*p00(x,prot) + sh*hl*p11(x,prot) + (1-sh)*hl*p01(x,prot) + sh*(1-hl)*p10(x,prot); }
    
  public:

    Curve (const FloatList& points);
   ~Curve ();

    float getVal (float x);
    void   getVal (const FloatList& t, FloatList& res);
};
};

#endif
